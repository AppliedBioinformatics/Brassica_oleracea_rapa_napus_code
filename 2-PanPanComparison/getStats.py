
import pandas as pd
from collections import OrderedDict

def check_napus():
    synths = set([x.rstrip() for x in open('../../synthetic_lines.txt')])

    main_dict = OrderedDict()# contains allelles for all napus inds
    main_dict_no_synths = OrderedDict() # contains alleles for all napus inds, without synths
    oler_dict = OrderedDict() # contains alleles for oleracea only
    rapa_dict = OrderedDict() # contains alleles for rapa only

    napus_core_out = open('Napus_Pangenome_With_Synths_Core_GeneNames.txt','w')
    napus_var_out = open('Napus_Pangenome_With_Synths_Variable_GeneNames.txt','w')

    napus_wo_synth_core_out = open('Napus_Pangenome_Without_Synths_Core_GeneNames.txt','w')
    napus_wo_synth_var_out = open('Napus_Pangenome_Without_Synths_Variable_GeneNames.txt','w')

    napus_only_napus_out = open('Napus_Pangenome_OnlyNapus.GeneNames.txt','w')

    with open('../1-FilterPavTables/NewPAV_Table_nice_names_intersection_of_individuals.filtered.csv') as fh:
        header = fh.readline().split()
        ind, genes = header[0], header[1:]
        print('Total genes: %s'%len(genes))
        for g in genes:
            main_dict[g] = []
            main_dict_no_synths[g] = []
            oler_dict[g] = []
            rapa_dict[g] = []
        core, variable = 0, 0
        core_without_synth, variable_without_synth = 0, 0
        rapa_present = 0
        oler_present = 0
        rapa_oler_present = 0
        neither_present = 0

        inds = []
        for line in fh:
            ll = line.split()
            napus = False
            if 'Brassica_napus' in ll[0]:
                napus = True
            oleracea = False
            if 'Brassica_oleracea' in ll[0]:
                oleracea = True
            rapa = False
            if 'Brassica_rapa' in ll[0]:
                rapa = True

            short_name = ll[0].replace('Brassica_napus_', '')
            syn = False
            if short_name in synths:
                syn = True

            for allele, g in zip(ll[1:], genes):
                if napus:
                    main_dict[g].append(allele)
                    # is this a filthy synth?
                    if not syn:
                        main_dict_no_synths[g].append(allele)
                if rapa:
                    rapa_dict[g].append(allele)
                if oleracea:
                    oler_dict[g].append(allele)

        for g in main_dict:
            alleles = main_dict[g]
            nosynth_alleles = main_dict_no_synths[g]

            variable_c = False
            for x in alleles:
                if x == '0':
                    variable_c = True
                    break

            if variable_c:
                variable += 1
                napus_var_out.write('%s\n'%g)
            else:
                core += 1
                napus_core_out.write('%s\n'%g)
                    
            # synths? 
            variable_s = False
            for x in nosynth_alleles:
                if x == '0':
                    variable_s = True
                    break
            if variable_s:
                variable_without_synth += 1
                napus_wo_synth_var_out.write('%s\n'%g)
            else:
                core_without_synth += 1
                napus_wo_synth_core_out.write('%s\n'%g)


            # now check across all species
            rapa_alleles = rapa_dict[g]
            oler_alleles = oler_dict[g]
            rapa_p = False
            oler_p = False
            oler_rapa_p = False
            neither = True
            for r, o in zip(rapa_alleles, oler_alleles):
                if r == '1' and o == '1':
                    oler_rapa_p = True
                if r == '1' and o == '0':
                    rapa_p = True
                if r == '0' and o == '1':
                    oler_p = True

            if oler_rapa_p:
                rapa_oler_present += 1
            elif rapa_p and not oler_p:
                rapa_present += 1
            elif oler_p and not rapa_p:
                oler_present += 1
            else:
                neither_present += 1
                napus_only_napus_out.write('%s\n'%g)


        print('Core, variable: %s, %s'%(core, variable))
        print('Without synths: Core, variable: %s, %s'%(core_without_synth, variable_without_synth))
        print('In oleracea and rapa: %s'%rapa_oler_present)
        print('In rapa: %s'%rapa_present)
        print('In oleracea: %s'%oler_present)
        print('Neither: %s'%neither_present)


# NOW COMES OLERACEA

def check_oler():
    main_dict = OrderedDict()# contains allelles for all napus inds
    main_dict_no_synths = OrderedDict() # contains alleles for all napus inds, without synths
    napus_dict = OrderedDict() # contains alleles for napus only
    rapa_dict = OrderedDict() # contains alleles for rapa only

    oleracea_core_out = open('Oleracea_Pangenome_Core_GeneNames.txt', 'w')
    oleracea_var_out = open('Oleracea_Pangenome_Variable_GeneNames.txt','w')

    oleracea_only_oleracea_out = open('Oleracea_Pangenome_OnlyOleracea.GeneNames.txt','w')
    with open('../1-FilterPAVTables/Oleracea_PAV_nice_names_intersection_of_individuals.filtered.csv') as fh:
        header = fh.readline().split()
        ind, genes = header[0], header[1:]
        print('Total genes: %s'%len(genes))
        for g in genes:
            main_dict[g] = []
            napus_dict[g] = []
            rapa_dict[g] = []
        core, variable = 0, 0
        rapa_present = 0
        napus_present = 0
        rapa_napus_present = 0
        neither_present = 0

        inds = []
        for line in fh:
            ll = line.split()
            napus = False
            if 'Brassica_napus' in ll[0]:
                napus = True
            oleracea = False
            if 'Brassica_oleracea' in ll[0]:
                oleracea = True
            rapa = False
            if 'Brassica_rapa' in ll[0]:
                rapa = True

            for allele, g in zip(ll[1:], genes):
                if oleracea:
                    main_dict[g].append(allele)
                if rapa:
                    rapa_dict[g].append(allele)
                if napus:
                    napus_dict[g].append(allele)

        for g in main_dict:
            alleles = main_dict[g]

            variable_c = False
            for x in alleles:
                if x == '0':
                    variable_c = True
                    break

            if variable_c:
                variable += 1
                oleracea_var_out.write('%s\n'%g)
            else:
                core += 1
                oleracea_core_out.write('%s\n'%g)
                    
            # now check across all species
            rapa_alleles = rapa_dict[g]
            napus_alleles = napus_dict[g]
            rapa_p = False
            napus_p = False
            napus_rapa_p = False
            neither = True
            for r, n in zip(rapa_alleles, napus_alleles):
                if r == '1' and n == '1':
                    napus_rapa_p = True
                if r == '1' and n == '0':
                    rapa_p = True
                if r == '0' and n == '1':
                    napus_p = True

            if napus_rapa_p:
                rapa_napus_present += 1
            elif rapa_p and not napus_p:
                rapa_present += 1
            elif napus_p and not rapa_p:
                napus_present += 1
            else:
                neither_present += 1
                oleracea_only_oleracea_out.write('%s\n'%g)


        print('Oleracea:')
        print('Core, variable: %s, %s'%(core, variable))
        print('In napus and rapa: %s'%rapa_napus_present)
        print('In rapa: %s'%rapa_present)
        print('In napus: %s'%napus_present)
        print('Neither: %s'%neither_present)







# now comes rapa
def check_rapa():
    main_dict = OrderedDict()# contains allelles for all napus inds
    main_dict_no_fpscs = OrderedDict() # contains alleles for all napus inds, without synths
    napus_dict = OrderedDict() # contains alleles for napus only
    oler_dict = OrderedDict() # contains alleles for rapa only

    fpscs = set([x.rstrip() for x in open('../../fpscs.txt')])

    rapa_no_fpscs_variables = open('Rapa_Pangenome_NoFPSc_Variable_GeneNames.txt','w')
    rapa_all_variables = open('Rapa_Pangenome_WithFPSc_Variable_GeneNames.txt','w')

    rapa_core_out = open('Rapa_Pangenome_WithFPSc_Core_GeneNames.txt','w')

    rapa_only_rapa_out = open('Rapa_Pangenome_OnlyRapa.GeneNames.txt','w')
    with open('../1-FilterPAVTables/Rapa_PAV_nice_names_intersection_of_individuals.filtered.csv') as fh:
        header = fh.readline().split()
        ind, genes = header[0], header[1:]
        print('Total genes: %s'%len(genes))
        for g in genes:
            main_dict[g] = []
            main_dict_no_fpscs[g] = []
            oler_dict[g] = []
            napus_dict[g] = []
        core, variable = 0, 0
        core_no_fpscs, variable_no_fpscs = 0, 0
        oler_present = 0
        napus_present = 0
        oler_napus_present = 0
        neither_present = 0

        inds = []
        for line in fh:
            ll = line.split()
            napus = False
            if 'Brassica_napus' in ll[0]:
                napus = True
            oleracea = False
            if 'Brassica_oleracea' in ll[0]:
                oleracea = True
            rapa = False
            fpsc = False
            if 'Brassica_rapa' in ll[0]:
                rapa = True
                short = ll[0].replace('Brassica_rapa_', '')
                if short in fpscs:
                    fpsc = True

            for allele, g in zip(ll[1:], genes):
                if rapa:
                    main_dict[g].append(allele)
                    # is this a filthy synth?
                    if not fpsc:
                        main_dict_no_fpscs[g].append(allele)
                if oleracea:
                    oler_dict[g].append(allele)
                if napus:
                    napus_dict[g].append(allele)

        for g in main_dict:
            alleles = main_dict[g]
            alleles_no_fpscs = main_dict_no_fpscs[g]

            variable_c = False
            for x in alleles:
                if x == '0':
                    variable_c = True
                    break

            if variable_c:
                variable += 1
                rapa_all_variables.write('%s\n'%g)
            else:
                core += 1
                rapa_core_out.write('%s\n'%g)
                    
            variable_f = False
            for x in alleles_no_fpscs:
                if x == '0':
                    variable_f = True
                    break
            if variable_f:
                variable_no_fpscs += 1
                rapa_no_fpscs_variables.write('%s\n'%g)
            else:
                core_no_fpscs += 1

            # now check across all species
            oler_alleles = oler_dict[g]
            napus_alleles = napus_dict[g]
            oler_p = False
            napus_p = False
            napus_oler_p = False
            neither = True
            for o, n in zip(oler_alleles, napus_alleles):
                if o == '1' and n == '1':
                    napus_oler_p = True
                if o == '1' and n == '0':
                    oler_p = True
                if o == '0' and n == '1':
                    napus_p = True

            if napus_oler_p:
                oler_napus_present += 1
            elif oler_p and not napus_p:
                oler_present += 1
            elif napus_p and not oler_p:
                napus_present += 1
            else:
                neither_present += 1
                rapa_only_rapa_out.write('%s\n'%g)


        print('Rapa:')
        print('Core, variable: %s, %s'%(core, variable))
        print('Core, variable no fpscs: %s, %s'%(core_no_fpscs, variable_no_fpscs))
        print('In napus and oleracea: %s'%oler_napus_present)
        print('In oleracea: %s'%oler_present)
        print('In napus: %s'%napus_present)
        print('Neither: %s'%neither_present)





check_napus()
check_oler()
check_rapa()
