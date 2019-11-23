from backend.models.synth import Synthesizer


class TestSynth:
    def test_given_an_adn_seq_synthesizes_the_seq_returning_the_correct_protein(self):
        albumin = 'ATGGATGCACACAAGAGTGAGGTTGCTCATCGGTTTAAAGATTTGGGAGAAGAAAATTTCAAAGCCTTGGTGTTGATTGCCTTTGCTCAGTATCTTC' \
                  'AGCAGTGTCCATTTGAAGATCATGTAAAATTAGTGAATGAAGTAACTGAATTTGCAAAAACATGTGTTGCTGATGAGTCAGCTGAAAATTGTGACAA' \
                  'ATCACTTCATACCCTTTTTGGAGACAAATTATGCACAGTTGCAACTCTTCGTGAAACCTATGGTGAAATGGCTGACTGCTGTGCAAAACAAGAACCT' \
                  'GAGAGAAATGAATGCTTCTTGCAACACAAAGATGACAATCCAAATCTCCCCCGATTGGTGAGACCAGAGGTTGATGTGATGTGCACTGCTTTTCATG' \
                  'ACAATGAAGAGACATTTTTGAAAAAATACTTATATGAAATTGCCAGAAGACATCCTTACTTTTATGCCCCGGAACTCCTTTTCTTTGCTAAAAGGTA' \
                  'TAAAGCTGCTTTTACAGAATGTTGCCAAGCTGCTGATAAAGCAGCCTGCCTGTTGCCAAAGCTCGATGAACTTCGGGATGAAGGGAAGGCTTCGTCT' \
                  'GCCAAACAGAGACTCAAGTGTGCCAGTCTCCAAAAATTTGGAGAAAGAGCTTTCAAAGCATGGGCAGTAGCTCGCCTGAGCCAGAGATTTCCCAAAG' \
                  'CTGAGTTTGCAGAAGTTTCCAAGTTAGTGACAGATCTTACCAAAGTCCACACGGAATGCTGCCATGGAGATCTGCTTGAATGTGCTGATGACAGGGC' \
                  'GGACCTTGCCAAGTATATCTGTGAAAATCAAGATTCGATCTCCAGTAAACTGAAGGAATGCTGTGAAAAACCTCTGTTGGAAAAATCCCACTGCATT' \
                  'GCCGAAGTGGAAAATGATGAGATGCCTGCTGACTTGCCTTCATTAGCGGCTGATTTTGTTGAAAGTAAGGATGTTTGCAAAAACTATGCTGAGGCAA' \
                  'AGGATGTCTTCTTGGGCATGTTTTTGTATGAATATGCAAGAAGGCATCCTGATTACTCTGTCGTACTGCTGCTGAGACTTGCCAAGACATATGAAAC' \
                  'CACTCTAGAGAAGTGCTGTGCCGCTGCAGATCCTCATGAATGCTATGCCAAAGTGTTCGATGAATTTAAACCTCTTATGGAAGAGCCTCAGAATTTA' \
                  'ATCAAACAAAATTGTGAGCTTTTTGAGCAGCTTGGAGAGTACAAATTCCAGAATGCGCTATTAGTTCGTTACACCAAGAAAGTACCCCAAGTGTCAA' \
                  'CTCCAACTCTTGTAGAGGTCTCAAGAAACCTAGGAAAAGTGGGCAGCAAATGTTGTAAACATCCTGAAGCAAAAAGAATGCCCTGTGCAGAAGACTA' \
                  'TCTATCCGTGGTCCTGAACCAGTTATGTGTGTTGCATGAGAAAACGCCAGTAAGTGACAGAGTCACCAAATGCTGCACAGAATCCTTGGTGAACAGG' \
                  'CGACCATGCTTTTCAGCTCTGGAAGTCGATGAAACATACGTTCCCAAAGAGTTTAATGCTGAAACATTCACCTTCCATGCAGATATATGCACACTTT' \
                  'CTGAGAAGGAGAGACAAATCAAGAAACAAACTGCACTTGTTGAGCTTGTGAAACACAAGCCCAAGGCAACAAAAGAGCAACTGAAAGCTGTTATGGA' \
                  'TGATTTCGCAGCTTTTGTAGAGAAGTGCTGCAAGGCTGACGATAAGGAAACCTGCTTTGCCGAGGAGGGTAAAAAACTTGTTGCTGCAAGTCAAGCT' \
                  'GCCTTAGGCTTATAACATCACATTTAAAAGCATCTCAGCCTACCATGAGAATAAGAGAAAGAAAATGAAGATCAAAAGCTTATTCATTCTGTTTTTC' \
                  'TTTTTCGTTGGTGTAAAAGCCAACACCCTGTCTAAAAAACATAAATTTCTTTAATCATTTTAATCATTTTGCCTCTTTTCTCTGTGCTTCAATTAAT' \
                  'AAAAAATGGAAAGAATCTAAAAAAACCCCCCCCCCCCCCCCCCTGCAGCAATAGCAACAACGTTGCGCAAACTATTAACTGGCGAA'

        seq = Synthesizer.accepting(Synthesizer.ADN, albumin).run()

        assert seq == 'MDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCA' \
                      'KQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDEL' \
                      'RDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKEC' \
                      'CEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYA' \
                      'KVFDEFKPLMEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLH' \
                      'EKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEK' \
                      'CCKADDKETCFAEEGKKLVAASQAALGL'

    def test_given_an_arn_seq_synthesizes_the_seq_returning_the_correct_protein(self):
        uvi1d = 'UUUUCCUAA'

        seq = Synthesizer.accepting(Synthesizer.ARN, uvi1d).run()

        assert seq == 'FS'

    def test_given_an_already_synthesized_protein_it_returns_the_same(self):
        protein = 'MDAHKSE'

        seq = Synthesizer.accepting(Synthesizer.PROTEIN, protein).run()

        assert seq == protein
