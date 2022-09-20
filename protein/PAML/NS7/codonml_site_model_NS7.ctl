      seqfile = /data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_45species_forPAML_longer.best.nuc.phy
     treefile = /data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_45sp_forPAML.txt
      outfile = /data/share/htp/TRNP1/paper_data/protein/PAML/NS7/NS7_model_output.txt

        noisy = 3 *how much rubbish on the screen(??)
      verbose = 1 *1:detailed output
      runmode = 0 *0:user defined tree; -2: pairwise, performs ML estimation of ds and dn in pairwise comparisons of protein                     *-coding sequences


      seqtype = 1 *1:codons
    CodonFreq = 2 *0:1/61 each for the standard genetic code; 1:from the av. nucl. freqs (F1X4); 
                  *2: from the average nucleotide frequencies at the three codon positions (F3X4); 3: free parameters

        ndata = 1
        clock = 0 *molecular clock states that nucleotide substitutions occur at a constant rate
                  *0:no clock, 1:clock; 2:local clock 
       aaDist = 0 *0: equal AA distances assumed; 1: Grantham's matrix



        model = 0 *2a This is the alternative model for the branch-site test of positive selection. The null model is
                  *also the branch-site model A but with ω2 = 1 fixed, specified by
                  *Model A1: model = 2, NSsites = 2, fix_omega = 1, omega = 1 



      NSsites = 7 *1: neutral selection; 2: selection; 7: beta, w ratio between 0,1; 8: beta+w, additional w


    fix_omega = 0

        getSE = 1

   Small_Diff = .5e-6
    cleandata = 0

       method = 0
