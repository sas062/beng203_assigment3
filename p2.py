import pandas as pd
import matplotlib.pyplot as plt

PROTON = 1.007276466812
WATER = 18.0153
CO = 27.99491462
NH3 = 17.02654910
CO2 = 43.98982924

# downloaded amino acid masses from https://proteomicsresource.washington.edu/protocols06/masses.php and condensed for simplicity
def load_amino_acid_masses(filename='amino_acid_masses.xlsx'):
    df = pd.read_excel(filename)
    return dict(zip(df['aa_code'], df['mass']))

def theoretical_mass(sequence, residue=False, aa_masses=None):
    # sequence is string of amino acids 
    # if residue is true calc residue mass instead of peptide mass
    if aa_masses is None:
        aa_masses = load_amino_acid_masses()

    # Calculate the mass of the peptide
    mass = sum(aa_masses[aa] for aa in sequence)

    # add water for full peptide mass
    if not residue:
        mass += WATER

    return mass

# part b calculation
def parse_mgf_precursor(filename):
    pepmass = None
    charge = None

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith('PEPMASS='):
                pepmass = float(line.split('=')[1].split()[0])
            elif line.startswith('CHARGE='):
                charge = int(line.split('=')[1].rstrip('+'))

            if pepmass is not None and charge is not None:
                break

    return pepmass, charge

def read_fasta(filename):
    sequence_lines = []

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('>'):
                sequence_lines.append(line)

    return ''.join(sequence_lines)

def ppm_error(observed_mass, theoretical_mass):
    return abs(observed_mass - theoretical_mass) / theoretical_mass * 1e6

# part c
def find_candidate_peptides(sequence, observed_mass, tol=30, aa_masses=None):
    prefix_masses = [0]
    for aa in sequence:
        prefix_masses.append(prefix_masses[-1] + aa_masses[aa])

    candidates = []
    for start in range(len(sequence)):
        for end in range(start + 1, len(sequence) + 1):
            peptide = sequence[start:end]
            theoretical = prefix_masses[end] - prefix_masses[start] + WATER
            error = ppm_error(observed_mass, theoretical)

            if error <= tol:
                candidates.append((peptide, start + 1, end, theoretical, error))

    return sorted(candidates, key=lambda candidate: candidate[-1])

# part d - assume monoisotopic masses
def ion_neutral_mass(fragment_residue_mass, ion_type):
    offsets = {
        'a': -CO,
        'b': 0,
        'c': NH3,
        'x': CO2,
        'y': WATER,
        'z': WATER - NH3,
    }

    return fragment_residue_mass + offsets[ion_type]

def generate_spectrum(peptide, aa_masses=None):
    if aa_masses is None:
        aa_masses = load_amino_acid_masses()
    
    n = len(peptide)
    ions = []

    for i in range(1,n):
        N = peptide[:i]
        C = peptide[i:]
        N_mass = sum(aa_masses[aa] for aa in N)
        C_mass = sum(aa_masses[aa] for aa in C)

        for ion_type in ['a', 'b', 'c']:
            mass = ion_neutral_mass(N_mass, ion_type)

            for charge in range(6):
                charge = charge + 1
                mz = (mass + charge * PROTON) / charge
                ions.append({
                    'ion_type': ion_type,
                    'ion_number': i,
                    'ion_label': f'{ion_type}{i}',
                    'fragment_sequence': N,
                    'charge': charge,
                    'mz': mz
                })
        
        c_ind = n-i

        for ion_type in ['x', 'y', 'z']:
            mass = ion_neutral_mass(C_mass, ion_type)

            for charge in range(6):
                charge = charge + 1
                mz = (mass + charge * PROTON) / charge
                ions.append({
                    'ion_type': ion_type,
                    'ion_number': c_ind,
                    'ion_label': f'{ion_type}{c_ind}',
                    'fragment_sequence': C,
                    'charge': charge,
                    'mz': mz
                })

    return ions

def parse_mgf_spectrum(filename):
    peaks = []

    with open(filename) as f:
        in_ions = False
        for line in f:
            line = line.strip()

            if line == "BEGIN IONS":
                in_ions = True
                continue
            elif line == "END IONS":
                break

            if in_ions and line and not line.startswith(("TITLE=", "RTINSECONDS=", "PEPMASS=", "CHARGE=")):
                parts = line.split()
                if len(parts) >= 2:
                    mz = float(parts[0])
                    intensity = float(parts[1])
                    peaks.append((mz, intensity))

    return peaks


def ppm_error_mz(observed_mz, theoretical_mz):
    return abs(observed_mz - theoretical_mz) / theoretical_mz * 1e6


def match_theoretical_to_experimental(theoretical_spectrum, experimental_peaks, tol_ppm=10):
    matches = []
    used_exp = set()

    for ion in theoretical_spectrum:
        best_idx = None
        best_peak = None
        best_error = float("inf")

        for idx, (exp_mz, exp_intensity) in enumerate(experimental_peaks):
            if idx in used_exp:
                continue

            err = ppm_error_mz(exp_mz, ion['mz'])
            if err <= tol_ppm and err < best_error:
                best_idx = idx
                best_peak = (exp_mz, exp_intensity)
                best_error = err

        if best_peak is not None:
            used_exp.add(best_idx)
            matches.append({
                'exp_mz': best_peak[0],
                'exp_intensity': best_peak[1],
                'theoretical_mz': ion['mz'],
                'ppm_error': best_error,
                'ion_type': ion['ion_type'],
                'ion_number': ion['ion_number'],
                'ion_label': ion['ion_label'],
                'fragment_sequence': ion['fragment_sequence'],
                'charge': ion['charge']
            })

    return matches


def plot_experimental_spectrum_with_matches(experimental_peaks, matches, peptide, output_file="annotated_spectrum.png"):
    exp_mz = [p[0] for p in experimental_peaks]
    exp_intensity_raw = [p[1] for p in experimental_peaks]

    max_intensity = max(exp_intensity_raw) if exp_intensity_raw else 1
    exp_intensity = [100 * i / max_intensity for i in exp_intensity_raw]

    top3 = sorted(matches, key=lambda x: x['exp_intensity'], reverse=True)[:3]

    top3_keys = {
        (m['exp_mz'], m['ion_label'], m['charge'], m['fragment_sequence'])
        for m in top3
    }

    other_matches = [
        m for m in matches
        if (m['exp_mz'], m['ion_label'], m['charge'], m['fragment_sequence']) not in top3_keys
    ]

    matched_mz = [m['exp_mz'] for m in other_matches]
    matched_intensity = [100 * m['exp_intensity'] / max_intensity for m in other_matches]

    top3_mz = [m['exp_mz'] for m in top3]
    top3_intensity = [100 * m['exp_intensity'] / max_intensity for m in top3]

    plt.figure(figsize=(16, 6))

    plt.vlines(exp_mz, [0], exp_intensity, color='lightgray', linewidth=1, label='Experimental peaks')

    if matched_mz:
        plt.vlines(matched_mz, [0], matched_intensity, color='red', linewidth=1.5, label='Matched peaks')

    if top3_mz:
        plt.vlines(top3_mz, [0], top3_intensity, color='blue', linewidth=2.5, label='Top 3 matched peaks')

    plt.xlabel("m/z")
    plt.ylabel("Relative Intensity")
    plt.title(f"Experimental Spectrum with Matched Ions\n{peptide}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    return top3

def main():
    aa_masses = load_amino_acid_masses()

    print(f"The theoretical residue mass is: {theoretical_mass('PEPTIDE', residue=True, aa_masses=aa_masses):.3f} Da")
    print(f"The theoretical peptide mass is: {theoretical_mass('PEPTIDE', aa_masses=aa_masses):.3f} Da")
    print(f"The ion mass for charge Z = 2 is: {(theoretical_mass('PEPTIDE', aa_masses=aa_masses) + (2 * PROTON)) / 2:.3f} Da")

    pepmass, charge = parse_mgf_precursor('CAH_test_01_scan1133.mgf')
    observed_peptide_mass = charge * pepmass - charge * PROTON
    observed_residue_mass = observed_peptide_mass - WATER

    print(f"\nPEPMASS: {pepmass}")
    print(f"Charge: {charge}+")
    print(f"Observed neutral peptide mass: {observed_peptide_mass:.4f} Da")
    print(f"Observed residue mass: {observed_residue_mass:.4f} Da")

    sequence = read_fasta('CAH_BOVIN_fasta.txt')
    candidates = find_candidate_peptides(sequence, observed_peptide_mass, tol=30, aa_masses=aa_masses)
    print("\nCandidate peptides within 30 ppm:")
    for peptide, start, end, theoretical, error in candidates:
        print(f"{peptide}\tpositions {start}-{end}\ttheoretical mass = {theoretical:.4f} Da\tppm error = {error:.2f}")

    candidates = find_candidate_peptides(sequence, observed_peptide_mass, tol=35, aa_masses=aa_masses)
    print("\nCandidate peptides within 35 ppm:")
    for peptide, start, end, theoretical, error in candidates:
        print(f"{peptide}\tpositions {start}-{end}\ttheoretical mass = {theoretical:.4f} Da\tppm error = {error:.2f}")

    for ind, candidate in enumerate(candidates, start=1):
        peptide, start, end, theoretical, error = candidate
        spectrum = generate_spectrum(peptide, aa_masses=aa_masses)

        aby_c2 = [ion for ion in spectrum if ion['ion_type'] in ['a', 'b', 'y'] and ion['charge'] == 2]
        ions_sorted = sorted(aby_c2, key=lambda ion: (ion['ion_type'], ion['ion_number']))

        df = pd.DataFrame(ions_sorted)
        df = df[['ion_label', 'mz']]
        df = df.rename(columns={'ion_label': 'Ion', 'mz': 'm/z'})
        df['m/z'] = df['m/z'].map(lambda x: f"{x:.4f}")
        df.to_csv(f"candidate_{ind}_spectrum.csv", index=False)

    best_peptide = candidates[0][0]

    experimental_peaks = parse_mgf_spectrum('CAH_test_01_scan1133.mgf')
    theoretical_spectrum = generate_spectrum(best_peptide, aa_masses=aa_masses)

    matches = match_theoretical_to_experimental(theoretical_spectrum, experimental_peaks)
    top3 = plot_experimental_spectrum_with_matches(experimental_peaks, matches, best_peptide, output_file="annotated_spectrum.png")

    print("\nTop 3 most intense matched peaks:")
    for match in top3:
        print(
            f"Experimental m/z = {match['exp_mz']:.4f}\t"
            f"Ion = {match['ion_label']}\t"
            f"Charge = {match['charge']}\t"
            f"Sequence = {match['fragment_sequence']}\t"
            f"PPM error = {match['ppm_error']:.2f}"
        )

if __name__ == "__main__":    
    main()
