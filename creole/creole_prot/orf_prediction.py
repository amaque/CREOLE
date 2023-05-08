#!/usr/bin/env python3.11
# encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import concurrent.futures
import math
from typing import Any

from . import helpers


class OrfsPrediction:
    """
    A class used to translate nucleotide to aminoacid

    Attributes:
        dic_codons (dict[str, str]): A codon table used to translate
            a genetic code into a sequence of amino acids
        codon_freq (dict[str, float]): Codon Usage Frequency Table(1)

    """

    dic_codons: dict[str, str] = helpers.dic_codons()
    codon_freq: dict[str, float] = helpers.codon_freq()

    def __init__(
        self,
        id_seq: str = 'read_001',
        sequence: str = 'ATCG',
        min_orf_len: int = 75,
    ) -> None:  # Minimal ORF length (nt)
        """
        The constructor for OrfsPrediction class

        Note:
            Do not include the `self` parameter in the ``Args`` section.

        Args:
            id_seq (str): The nucleotide sequence header of the fasta or
                fastq file
            sequence (str): The nucleotide sequence
            min_orf_len (int): The minimum acceptable ORF length (by default 75
                nucleotide)
        """
        sequence = (
            sequence.upper()
            .replace('U', 'T')
            .replace('\n', '')
            .replace('N', '')
        )
        self.sequence: str = sequence
        self.sequenceLen: int = len(self.sequence)
        self.idSequence: str = id_seq
        self.seqType: str = 'nucleotide'
        self.minORFsLen: int = min_orf_len
        self.valid_sequence: dict[str, list[str]] = {
            'nucleotide': ['A', 'T', 'C', 'G']
        }
        self.is_valid: bool = self.__validation()
        assert self.is_valid, (
            f'Provided data does not seem to be a correct '
            f'{self.seqType} sequence'
        )

    def __validation(self) -> bool:
        """Check if it's a valid sequences"""
        return set(self.valid_sequence[self.seqType]).issuperset(self.sequence)

    def reverse_complement(self) -> str:
        """
        Swapping adenine with thymine and guanine with cytosine.
        Reversing newly generated string
        """
        remapping = str.maketrans('ATCG', 'TAGC')
        return self.sequence.translate(remapping)[::-1]

    def orfs_translation(
        self, start_pos: int = 0
    ) -> tuple[list[str], list[Any], list[int], list[str]]:
        """Translates nucleotide sequence into an aminoacid sequence"""
        amino_acid = [
            OrfsPrediction.dic_codons[self.sequence[pos : pos + 3]]
            for pos in range(start_pos, len(self.sequence) - 2, 3)
        ]
        freq = [
            OrfsPrediction.codon_freq[self.sequence[pos : pos + 3]]
            for pos in range(start_pos, len(self.sequence) - 2, 3)
        ]
        coord = [
            pos + 3 for pos in range(start_pos, len(self.sequence) - 2, 3)
        ]
        nt = [
            self.sequence[pos : pos + 3]
            for pos in range(start_pos, len(self.sequence) - 2, 3)
        ]
        return amino_acid, freq, coord, nt

    def nucleotide2aminoacid(self, start_pos: int = 0) -> list[str]:
        """Translates nucleotide sequence into an aminoacid sequence"""
        amino_acid = [
            OrfsPrediction.dic_codons[self.sequence[pos : pos + 3]]
            for pos in range(start_pos, len(self.sequence) - 2, 3)
        ]
        return amino_acid

    def generate_orfs(
        self, strand: str = 'direct', one: bool = False
    ) -> tuple[list[str], list[str], list[str]] | Any:
        """
        strand="both" (Default) :  Generate the six reading frames of a
        nucleotide sequence, including reverse complement
        strand="direct" : Generate the three reading frames of a
        nucleotide sequence, forward (+)
        strand="reverse" : Generate the three reading frames of a
        nucleotide sequence, reverse (-)
        """
        reverse = OrfsPrediction(sequence=self.reverse_complement())
        if strand == 'direct':
            frames1, frames2, frames3 = [], [], []
            with concurrent.futures.ThreadPoolExecutor() as executor:
                executor.submit(frames1.append, self.nucleotide2aminoacid(0))
                executor.submit(frames2.append, self.nucleotide2aminoacid(1))
                executor.submit(frames3.append, self.nucleotide2aminoacid(2))
            return frames1, frames2, frames3
        elif strand == 'both':
            forward_f1, forward_f2, forward_f3 = [], [], []
            reverse_f1, reverse_f2, reverse_f3 = [], [], []
            with concurrent.futures.ThreadPoolExecutor() as executor:
                executor.submit(
                    forward_f1.append, self.nucleotide2aminoacid(0)
                )
                executor.submit(
                    forward_f2.append, self.nucleotide2aminoacid(1)
                )
                executor.submit(
                    forward_f3.append, self.nucleotide2aminoacid(2)
                )
                executor.submit(
                    reverse_f1.append, reverse.nucleotide2aminoacid(0)
                )
                executor.submit(
                    reverse_f2.append, reverse.nucleotide2aminoacid(1)
                )
                executor.submit(
                    reverse_f3.append, reverse.nucleotide2aminoacid(2)
                )
            return (forward_f1, forward_f2, forward_f3), (
                reverse_f1,
                reverse_f2,
                reverse_f3,
            )
        elif strand == 'reverse':  # Validate the strand in the main function
            frames1 = []
            frames2 = []
            frames3 = []
            with concurrent.futures.ThreadPoolExecutor() as executor:
                executor.submit(
                    frames1.append, reverse.nucleotide2aminoacid(0)
                )
                executor.submit(
                    frames2.append, reverse.nucleotide2aminoacid(1)
                )
                executor.submit(
                    frames3.append, reverse.nucleotide2aminoacid(2)
                )
            return frames1, frames2, frames3
        del reverse

    def possible_proteins_def1(self, aminoacid_seqs):
        """
        Compute all possible proteins in an aminoacid seq and return a list of
        possible proteins
        Definition 1: an ORF is a sequence that has a length divisible by three
        and begins with a translation start
        codon (ATG) and ends at a stop codon
        """
        protein = []
        proteins = []
        coord = []
        coords = []
        freq = []
        freqs = []
        prob = 1 / 64
        n = []
        for aminoacid_seq in aminoacid_seqs:
            for aminoacid, codon_frq, coords_ in zip(
                aminoacid_seq[0], aminoacid_seq[1], aminoacid_seq[2]
            ):
                # print("@@",coords_)
                if aminoacid == '_':
                    # STOP accumulating amino acids if _ - STOP was found
                    if protein:
                        for p, f, c in zip(protein, freq, coord):
                            if p:
                                proteins.append(p)
                                # if len(p) <= 1193:
                                freqs.append(f)  # - (prob ** len(p)))
                                coords.append((c, c + (len(p) * 3) - 3))
                        protein = []
                        freq = []
                        coord = []
                else:
                    # START accumulating amino acids if M - START was found
                    if aminoacid == 'M':
                        protein.append('')
                        freq.append(0)
                        coord.append(coords_)
                    for i in range(len(protein)):
                        protein[i] += aminoacid
                        freq[i] += math.log10(codon_frq / prob)
                        # coord.append(coords_)
        if proteins:
            return proteins, freqs, coords

    def possible_proteins_def2(self, aminoacid_seqs):
        """
        Compute all possible proteins in an aminoacid seq and return a list of
        possible proteins
        Definition 2: an ORF is a sequence that has a length divisible by three
        and is bounded by stop codons
        """
        protein = ''
        proteins = list()
        freq = 1
        freqs = list()
        prob = 1 / 64
        coord = []
        add_coord = list()

        for aminoacid_seq in aminoacid_seqs:
            for aminoacid, codon_frq, coords in zip(
                aminoacid_seq[0], aminoacid_seq[1], aminoacid_seq[2]
            ):
                if aminoacid != '_' and len(protein) <= len(aminoacid_seq[0]):
                    protein += aminoacid
                    # print("Frequency: ", freq)

                    freq += math.log10(codon_frq / prob)
                    add_coord.append(coords)
                else:
                    if protein:
                        try:
                            if len(protein) >= self.minORFsLen / 3:
                                proteins.append(protein)
                                freqs.append(
                                    freq - 1
                                )  # - (prob ** len(protein)))
                                coord.append((add_coord[0], add_coord[-1]))
                        except IndexError as err:
                            print(protein, err)
                    protein = ''
                    freq = 0
                    add_coord = []
        return proteins, freqs, coord

    def possible_proteins_def2_backup(self, aminoacid_seq):
        """
        Compute all possible proteins in an aminoacid seq and return a list of
        possible proteins
        Definition 2: an ORF is a sequence that has a length divisible by three
         and is bounded by stop codons
        """
        protein = ''
        proteins = []

        for aminoacid in aminoacid_seq:
            if aminoacid != '_' and len(protein) <= len(aminoacid_seq):
                protein += aminoacid
            else:
                if len(protein) >= self.minORFsLen / 3:
                    proteins.append(protein)
                protein = ''
        return proteins

    def all_proteins_from_orfs(
        self, start_pos=0, end_pos=0, ordered=False, one=False
    ):
        """Compute all possible proteins for all open reading frames"""
        if one:
            if end_pos > start_pos:
                tmp_seq = OrfsPrediction(
                    self.sequence[start_pos:end_pos], self.minORFsLen
                )
                orfs = tmp_seq.generate_orfs(one=True)
            else:
                orfs = self.generate_orfs(one=True)

        else:
            if end_pos > start_pos:
                tmp_seq = OrfsPrediction(
                    self.sequence[start_pos:end_pos], self.minORFsLen
                )
                orfs = tmp_seq.generate_orfs()
            else:
                orfs = self.generate_orfs()
        # print("-----------",orfs)
        protein_list = []
        for orf in orfs[0]:
            with concurrent.futures.ThreadPoolExecutor() as executor:
                # def1 = executor.submit(self.possible_proteins_def1, orf)
                def2 = executor.submit(self.possible_proteins_def2, orf)
            proteins = def2.result()  # + def2.result()

            for protein in proteins:
                if protein:
                    protein_list.append(protein)
        # reverse
        reverse = []
        for orf in orfs[1]:
            with concurrent.futures.ThreadPoolExecutor() as executor:
                # def1 = executor.submit(self.possible_proteins_def1, orf)
                def2 = executor.submit(self.possible_proteins_def2, orf)
            reverse_proteins = def2.result()  # + def2.result()

            for protein in reverse_proteins:
                if protein:
                    # print (protein)
                    reverse.append(protein)

        # a = sorted(protein_list, key=len, reverse=True)
        # a = protein_list
        # print("BBBBBBBBBBBBBBBBB ",protein_list)
        # print("KKKKKKKKKKKKKKKKK ",a)

        # count = 1
        # for k in a:
        #    print("Frame ",count,":", k)
        #    count += 1

        if ordered:
            return sorted(protein_list, key=len, reverse=True)
        return protein_list, reverse

    def all_proteins_from_orfs_1(self, start_pos=0, end_pos=0, ordered=False):
        """Compute all possible proteins for all open reading frames"""
        if end_pos > start_pos:
            tmp_seq = OrfsPrediction(
                self.sequence[start_pos:end_pos], self.minORFsLen
            )
            orfs = tmp_seq.generate_orfs()
        else:
            orfs = self.generate_orfs()

        protein_list = []
        for orf in orfs:
            proteins = self.possible_proteins_def2(orf)
            for protein in proteins:
                protein_list.append(protein)

        if ordered:
            return sorted(protein_list, key=len, reverse=True)
        return protein_list

    def longer_protein(self):
        """Scroll through all possible proteins and select the longest one"""
        orfs = self.all_proteins_from_orfs()
        if orfs:
            check = False
            if str(orfs[0]).startswith('M'):
                return orfs[0]
            else:
                for orf in orfs:
                    if orf.startswith('M') and len(orfs[0]) - len(orf) <= 3:
                        check = True
                        return orf
            if str(orfs[0])[0] != 'M' and check is False:
                return orfs[0]

    def multi_longer_protein(self):
        """Scroll through all possible proteins and select the longest ones"""
        if self.all_proteins_from_orfs():
            multi_orf = list()
            orf_list = sorted(
                self.all_proteins_from_orfs(), key=len, reverse=True
            )
            for i in range(len(orf_list)):
                if len(orf_list[0]) == len(orf_list[i]):
                    multi_orf.append(orf_list[i])
                else:
                    break
            return multi_orf

    def best_orf(self):
        """select the best ORFs"""
        protein = []
        score = []
        coord = list()
        frame = []
        orf_config = dict()
        g_dic = dict()
        count = 0

        if self.all_proteins_from_orfs():
            for prots in self.all_proteins_from_orfs():
                if prots:
                    count = count + 1
                    for prot in prots:
                        if type(prot) == str:
                            protein.append(prot)
                            frame.append(count)
                        elif type(prot) == tuple:
                            coord.append(prot)
                        elif isinstance(prot, float):
                            score.append(prot)
                        # elif  isinstance(prot, int):
                        #     frame.append(count)

            # print("Prot ", protein)
            # print("coord ", coord)
            # print("score ", score)
            # print("Frame ", frame)

        # print("Prot ", g_dic.index())
        # print("coord ", g_dic.value()[2])
        # print("score ", g_dic.value()[2])

        for p, c, s, f in zip(protein, coord, score, frame):
            orf_config[p] = [c, s, len(p), f]

        # print("SORTED",sorted(protein,key=len, reverse=True))
        # print(coord[protein.index("TRPSGGDIHTVESLRFHTEKNHSPFSEAQDTPNRYHRVHGTSKC")])

        if score:
            best_protein = protein[score.index(max(score))]
            if best_protein is not None and max(score) > 0:
                finalOrfs = (
                    self.idSequence
                    + ' sc='
                    + str(max(score))
                    + ' len='
                    + str(len(best_protein)),
                    best_protein,
                )
                return finalOrfs
            else:
                print(self.idSequence, 'no protein')

    def orfs_list(self, reverse=0):
        """ORFs info: Aminoacids, Frames, Coordinates"""
        protein = []
        score = []
        coord = list()
        count = 0
        frame = []
        orfs = self.all_proteins_from_orfs()[reverse]
        if orfs:
            for prots in orfs:
                if prots:
                    count = count + 1
                    for prot in prots:
                        if type(prot) == str:
                            protein.append(prot)
                            frame.append(count)
                        elif type(prot) == tuple:
                            coord.append(prot)
                        elif isinstance(prot, float):
                            score.append(prot)
            i = 0
            while i < len(frame):
                if frame[i] == 4:
                    frame[i] = 2
                if frame[i] == 7:
                    frame[i] = 3
                if frame[i] == 10:
                    frame[i] = 1
                if frame[i] == 13:
                    frame[i] = 2
                if frame[i] == 16:
                    frame[i] = 3
                i += 1
        if protein:
            if reverse == 0:
                # print("Teste retorno",protein, coord, score, frame, self.idSequence, self.sequence)
                return (
                    protein,
                    coord,
                    score,
                    frame,
                    self.idSequence,
                    self.sequence,
                )
            else:
                return (
                    protein,
                    coord,
                    score,
                    frame,
                    self.idSequence,
                    self.reverse_complement(),
                )

    def run_orfs_list(self, reverse=0):
        if reverse == 0:
            return self.orfs_list(reverse=0)  # , self.orfs_list(reverse=1)
        elif reverse == 1:
            return self.orfs_list(reverse=0), self.orfs_list(reverse=1)

    def assemble_orf(self):
        """gather the best ORFs"""
        protein = []
        score = []
        coord = list()
        orf_config = dict()
        count = 0
        frame = []
        if self.all_proteins_from_orfs():
            for prots in self.all_proteins_from_orfs():
                if prots:
                    count = count + 1
                    for prot in prots:
                        if type(prot) == str:
                            protein.append(prot)
                            frame.append(count)
                        elif type(prot) == tuple:
                            coord.append(prot)
                        elif isinstance(prot, float):
                            score.append(prot)

            print('Prot ', protein)
            print('coord ', coord)
            print('score ', score)
            print('Frame ', frame)
            i = 0
            while i < len(frame):
                # replace 4 with 2
                if frame[i] == 4:
                    frame[i] = 2
                # replace 7 with 3
                if frame[i] == 7:
                    frame[i] = 3
                i += 1
            print('Frame ', frame)

        protein_sort = sorted(protein, key=len, reverse=True)
        assemb_prot = []
        start = 0
        end = 0
        count = 1

        def insert_before(base, new):
            compliment_start = coord[protein.index(new)][0]
            base_start = coord[protein.index(base)][0]
            compliment_end = coord[protein.index(new)][1]
            base_end = coord[protein.index(base)][1]

            print(compliment_start, compliment_end, base_start, base_end)
            if compliment_start < base_start:
                if compliment_end < base_end:
                    if (
                        compliment_end - base_start
                    ) / 2 + compliment_start > 44:
                        if (compliment_end - base_start) / 2 + base_end > 44:
                            print(int((compliment_end - base_start) / 4))
                            basex = base[0][
                                int((compliment_end - base_start) / 4) :
                            ]
                            newx = new[0][
                                : -int((compliment_end - base_start) / 4)
                            ]
                            print('basex: ', basex)
                            print('Newx: ', newx)
                            return newx + basex
                        else:
                            print('rule 4')
                    else:
                        print('rule 3')
                else:
                    print('rule 2')
            else:
                print('rule 1')

        def insert_after(base, compliment):
            print('PPPPPPPPPPPPP: ', base)

            if (
                coord[protein.index(base)][0]
                > coord[protein.index(compliment)][0]
            ):
                a = base
                b = compliment
                base = b
                compliment = a

            compliment_start = coord[protein.index(compliment)][0]
            base_start = coord[protein.index(base)][0]
            compliment_end = coord[protein.index(compliment)][1]
            base_end = coord[protein.index(base)][1]
            compliment_frame = frame[protein.index(compliment)]
            base_frame = frame[protein.index(base)]

            base_start_inter = 0
            compliment_end_inter = 0
            if base_frame == 1 and compliment_frame == 2:
                base_start_inter = compliment_start - 1
            elif base_frame == 2 and compliment_frame == 1:
                base_start_inter = compliment_start - 2
                compliment_end_inter = base_end - 2
            elif base_frame == 2 and compliment_frame == 3:
                base_start_inter = compliment_start
                compliment_end_inter = base_end - 3

            print('FRAME: ', base_frame, compliment_frame)
            print('JUMB: ', base_start_inter)

            print(compliment_start, compliment_end, base_start, base_end)

            if compliment_start > base_start:
                if compliment_end > base_end:
                    if (
                        len(compliment) > 44
                        and compliment_start - base_start > 44
                    ):
                        if compliment_end - base_end > 44:
                            print(int((base_end - compliment_start) / 2))

                            interc = int((base_end - compliment_start) / 2)

                            # basex =[OrfsPrediction.dic_codons[self.sequence[pos:pos + 3]] for pos in
                            #               range(base_start - 3, len(self.sequence[:base_end-interc]) - 2, 3)]

                            baseSP = [
                                OrfsPrediction.dic_codons[
                                    self.sequence[pos : pos + 3]
                                ]
                                for pos in range(
                                    base_start_inter,
                                    len(self.sequence[:base_end]) - 2,
                                    3,
                                )
                            ]

                            complimentSP = [
                                OrfsPrediction.dic_codons[
                                    self.sequence[pos : pos + 3]
                                ]
                                for pos in range(
                                    compliment_start,
                                    len(self.sequence[:base_end]) - 2,
                                    3,
                                )
                            ]

                            print(self.sequence[base_start_inter:base_end])
                            print(self.sequence[compliment_start:base_end])
                            # for x in self.sequence[base_start_inter:base_end + 3]:
                            dic_f = {}
                            for x in range(
                                base_start_inter, len(self.sequence[:base_end])
                            ):

                                # print (self.sequence[base_start_inter:x])

                                xx = self.sequence[base_start_inter - 1 : x]
                                xy = self.sequence[x - 2 : base_end - 4]
                                inter = xx + xy

                                # print("INTERRRRRR: ", xy)

                                frq = [
                                    OrfsPrediction.codon_freq[
                                        inter[pos : pos + 3]
                                    ]
                                    for pos in range(0, len(inter[:x]) - 2, 3)
                                ]
                                frequency = 0
                                prob = 1 / 64
                                FQ = []

                                for f in frq:
                                    frequency += math.log10(f / prob)
                                FQ.append(frequency)

                                dic_f[
                                    ''.join(
                                        [
                                            OrfsPrediction.dic_codons[
                                                inter[pos : pos + 3]
                                            ]
                                            for pos in range(
                                                0, len(inter[:x]) - 2, 3
                                            )
                                        ]
                                    )
                                ] = [a for a in FQ][0]

                                # print("".join([OrfsPrediction.dic_codons[inter[pos:pos + 3]] for pos in
                                #               range(0, len(inter[:x]) - 2, 3)]), [a for a in FQ])

                            for x, y in dic_f.items():
                                print(x, y)

                            print(
                                'maior: ',
                                max(dic_f, key=dic_f.get),
                                max(dic_f.values()),
                            )
                            itc = max(dic_f, key=dic_f.get)
                            print(
                                base[: -len(itc)]
                                + '-'
                                + itc
                                + '-'
                                + compliment[len(itc) :]
                            )

                            # print([OrfsPrediction.dic_codons[self.sequence[pos:pos + 3]] for pos in
                            #               range(compliment_start, len(self.sequence[x:]) - 2, 3)])
                            # print("XXXXX: ", xx)

                            # print(xx + xy)

                            # v = OrfsPrediction(idSeq="01", sequence= xx + xy, minorfslen=3)
                            # vv = v.all_proteins_from_orfs(one=True)
                            # print(v.best_orf())
                            # print(vv)
                            # print(xx, xy)
                            # print("XXXXX: ", xy)
                            # k = x + self.sequence[compliment_start + 3:base_end]
                            # print(k)

                            # o = [self.sequence[pos:pos + 3] for pos in
                            #               range(base_start_inter, len(self.sequence[:base_end]) - 2, 3)]
                            #
                            # n = [self.sequence[pos:pos + 3] for pos in
                            #               range(compliment_start, len(self.sequence[:base_end]) - 2, 3)]

                            # print(o)
                            # print(n)
                            # complimentx =[OrfsPrediction.dic_codons[self.sequence[pos:pos + 3]] for pos in
                            #               range((compliment_start+4) + interc, len(self.sequence[:compliment_end]) - 2, 3)]

                            # basex = base[0][:- int((base_end - compliment_start)/2)]

                            # complimentx = compliment[0][int((base_end - compliment_start)/2):]
                            # print("basex: ", basex)
                            # print("complimentx: ", complimentx)
                            print('baseSP: ', baseSP, len(baseSP))
                            print(
                                'complimentSP: ',
                                complimentSP,
                                len(complimentSP),
                            )
                            # return basex+complimentx
                        else:
                            print('rule 4')
                    else:
                        print('rule 3')
                else:
                    print('rule 2')
            else:
                print('rule 1')

        print('Mais longo', protein_sort[0])
        print('2 Mais longo', protein_sort[1])
        print(insert_after(protein_sort[0], protein_sort[1]))

        for p in protein_sort:
            print(
                'ORF%d:' % count,
                p,
                coord[protein.index(p)][0],
                coord[protein.index(p)][1],
            )  # , score[protein.index(p)])
            if assemb_prot:
                if coord[protein.index(p)][1] > end and len(p) > 40:
                    assemb_prot.append(p)
                    end = coord[protein.index(p)][1]
                elif (
                    coord[protein.index(p)][1] <= start
                    and start - coord[protein.index(p)][0] > 0
                ):
                    assemb_prot.insert(0, p)
                    start = coord[protein.index(p)][0]
            else:
                assemb_prot.append(p)
                end = coord[protein.index(p)][1]
                start = coord[protein.index(p)][0]
            count += 1

        print('assemb_prot:', assemb_prot, self.sequenceLen / 3)

        # for p,c,s in zip(protein,coord,score):
        #
        #     orf_config[p] = [c,s, len(p)]
        #
        # print(orf_config)
        #
        # if score:
        #     best_protein = protein[score.index(max(score))]
        #     if best_protein is not None and max(score) > 0:
        #         finalOrfs = (
        #             self.idSequence + " sc=" + str(max(score)) + " len=" + str(len(best_protein)), best_protein)
        #         return finalOrfs
        #     else:
        #         print(self.idSequence, "no protein")

    def fasta2dicionary(self, file_path):
        """Read fasta"""
        with open(file_path, 'r') as rfile:
            fasta = [l.strip() for l in rfile.readlines()]

        dic_fasta = {}
        ids = ''

        for line in fasta:
            if '>' in line:
                ids = line
                dic_fasta[ids] = ''
            else:
                dic_fasta[ids] += line
        return dic_fasta
