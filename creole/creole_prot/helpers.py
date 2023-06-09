def dic_codons() -> dict[str, str]:
    """
    human codon
    :return: dict string, string
    """
    return dict(
        GCT='A',
        GCC='A',
        GCA='A',
        GCG='A',
        TGT='C',
        TGC='C',
        GAT='D',
        GAC='D',
        GAA='E',
        GAG='E',
        TTT='F',
        TTC='F',
        GGT='G',
        GGC='G',
        GGA='G',
        GGG='G',
        CAT='H',
        CAC='H',
        ATA='I',
        ATT='I',
        ATC='I',
        AAA='K',
        AAG='K',
        TTA='L',
        TTG='L',
        CTT='L',
        CTC='L',
        CTA='L',
        CTG='L',
        ATG='M',
        AAT='N',
        AAC='N',
        CCT='P',
        CCC='P',
        CCA='P',
        CCG='P',
        CAA='Q',
        CAG='Q',
        CGT='R',
        CGC='R',
        CGA='R',
        CGG='R',
        AGA='R',
        AGG='R',
        TCT='S',
        TCC='S',
        TCA='S',
        TCG='S',
        AGT='S',
        AGC='S',
        ACT='T',
        ACC='T',
        ACA='T',
        ACG='T',
        GTT='V',
        GTC='V',
        GTA='V',
        GTG='V',
        TGG='W',
        TAT='Y',
        TAC='Y',
        TAA='_',
        TAG='_',
        TGA='_',
    )


def codon_freq() -> dict[str, float]:
    """
    human codon frequency
    :return: dict string, float
    """
    return dict(
        TTT=0.46,
        TCT=0.19,
        TAT=0.44,
        TGT=0.46,
        TTC=0.54,
        TCC=0.22,
        TAC=0.56,
        TGC=0.54,
        TTA=0.08,
        TCA=0.15,
        TAA=0.30,
        TGA=0.47,
        TTG=0.13,
        TCG=0.05,
        TAG=0.24,
        TGG=0.99,
        CTT=0.13,
        CCT=0.29,
        CAT=0.42,
        CGT=0.08,
        CTC=0.20,
        CCC=0.32,
        CAC=0.58,
        CGC=0.18,
        CTA=0.07,
        CCA=0.28,
        CAA=0.27,
        CGA=0.11,
        CTG=0.40,
        CCG=0.11,
        CAG=0.73,
        CGG=0.20,
        ATT=0.36,
        ACT=0.25,
        AAT=0.47,
        AGT=0.15,
        ATC=0.47,
        ACC=0.36,
        AAC=0.53,
        AGC=0.24,
        ATA=0.17,
        ACA=0.28,
        AAA=0.43,
        AGA=0.21,
        ATG=0.99,
        ACG=0.11,
        AAG=0.57,
        AGG=0.21,
        GTT=0.18,
        GCT=0.27,
        GAT=0.46,
        GGT=0.16,
        GTC=0.24,
        GCC=0.40,
        GAC=0.54,
        GGC=0.34,
        GTA=0.12,
        GCA=0.23,
        GAA=0.42,
        GGA=0.25,
        GTG=0.46,
        GCG=0.11,
        GAG=0.58,
        GGG=0.25,
    )
