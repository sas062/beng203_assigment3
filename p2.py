import pandas as pd

def theoretical_mass(sequence, residue=False):
    # sequence is string of amino acids 
    # if residue is true calc residue mass instead of peptide mass
    df = pd.read_excel('amino_acid_masses.xlsx')
    aa_masses = dict(zip(df['aa_code'], df['mass']))

    # Calculate the mass of the peptide
    mass = sum(aa_masses[aa] for aa in sequence)

    # add water for full peptide mass
    if not residue:
        mass += 18.0153

    return mass

def main():
    print(f"The theoretical residue mass is: {theoretical_mass('PEPTIDE', residue=True):.3f} Da")
    print(f"The theoretical peptide mass is: {theoretical_mass('PEPTIDE'):.3f} Da")
    print(f"The ion mass for charge Z = 2 is: {(theoretical_mass('PEPTIDE') + 2) / 2:.3f} Da")


if __name__ == "__main__":    
    main()