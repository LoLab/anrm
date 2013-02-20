"""
    Overview
    ========
    
    PySB implementations of the apoptosis-necrosis reaction model version 1.0
    (ANRM 1.0) originally published in [Irvin,NECRO2013]_.
    
    This file contains functions that implement the extrinsic apoptosis pathway
    in five modules:
    
    - CD95 Ligation to formation of secondary complex
    - TNFR1 ligation to formation of complex II
    - Secondary complexes to formation of Riptosomes and Necrosomes and Bid
    activation.
    - Execution of apoptosis (MOMP from Albeck)
    - Execution of necrosis
    
    For the (MOMP) segment there are five variants, which correspond to the five
    models described in Figure 11 of [Albeck2008]_:
    
    - "Minimal Model" (Figure 11b, :py:func:`albeck_11b`)
    - "Model B + Bax multimerization" (Figure 11c, :py:func:`albeck_11c`)
    - "Model C + mitochondrial transport" (Figure 11d, :py:func:`albeck_11d`)
    - "Current model" (Figure 11e, :py:func:`albeck_11e`)
    - "Current model + cooperativity" (Figure 11f, :py:func:`albeck_11f`)
    """

import numpy
from pysb import *
from pysb.util import alias_model_components
from pysb.macros import *

#from shared_anrm import *
#from earm.shared import *

Model()
    
Parameter('KF', 1e-6) # Generic association rate constant
Parameter('KR', 1e-3) # Generic dessociation rate constant
Parameter('KC', 1)    # Generic catalytic rate constant
Parameter('KC2', 10) # Generic catalytic rate constant
Parameter('KE', 1e-4) # Generic gene expression rate constant

Parameter('Ka_RIP1_FADD',   1e-7) # Biochemica et Biophysica Acta 1834(2013) 292-300
Parameter('Kd_RIP1_FADD',   1e-8) # Biochemica et Biophysica Acta 1834(2013) 292-300

Parameter('Kf_Apaf_acti',   5e-7) # from Albeck_modules.py
Parameter('Kf_Apop_asse',   5e-8) # from Albeck_modules.py
Parameter('Kf_C3_activa',   5e-9) # from Albeck_modules.py
Parameter('Kf_Apop_inhi',   2e-6) # from Albeck_modules.py
Parameter('Kf_Smac_inhi',   7e-6) # from Albeck_modules.py
Parameter('Kf_C3_activ2',   1e-7) # from Albeck_modules.py
Parameter('Kf_C3_ubiqui',   2e-6) # from Albeck_modules.py
Parameter('Kc_C3_ubiqui',   1e-1) # from Albeck_modules.py
Parameter('Kr_PARP_clea',   1e-2) # from Albeck_modules.py
Parameter('Kf_C8_activ2',   3e-8) # from Albeck_modules.py
Parameter('Kf_Bax_activ',   1e-7) # from Albeck_modules.py

def CD95_to_SecondaryComplex_monomers():
    """ Declares Fas ligand, CD95, FADD, Flip_L, Flip_S procaspase8 and Caspase 8.
        
    'bf' is the site to be used for all binding reactions.
        
    The 'state' site denotes various localization and/or activity states of a
    Monomer, with 'C' denoting cytoplasmic localization and 'M' mitochondrial
    localization.
    """
    Monomer('Fas', ['blig'])    #Fas Ligand
    Monomer('CD95', ['blig', 'bDD'])           #Fas Receptor (CD95)
    Monomer('FADD', ['bDD', 'bDED1','bDED2'])    #FADD
    Monomer('flip_L', ['bDED'])   #c-Flip[L] binds FADD at bca1 or bca2
    Monomer('flip_S', ['bDED'])   #c-Flip[S] binds FADD at bca1 or bca2
    Monomer('proC8', ['bDED'])    #procaspase 8 binds FADD at bca1 or bca2
    Monomer('C8', ['bC8'])        #active caspase 8
    Monomer('Bar', ['bC8'])       #bifunctional apoptosis regulator

def CD95_to_SecondaryComplex():
    """Defines the interactoins from CD95 ligation to generation of secondary
    complexes as per ANRM 1.0.
    
    Uses Fas, CD95, FADD, flip_L, flip_S, proC8 monomers and C8 active dimers, and
    their assicated parameters to generate rules that describe the ligand/receptor
    binding, FADD recruitment, proC8 and c-flip recruitment, activation of caspase
    and release of the secondary complex. This model assumes that one copy of proC8
    binds FADD before c-flip. 
    
    This model converts proC8:proC8 to C8 (active caspase 8 dimer)
    This model also produces Secondary complex, FADD:proC8:c-Flip.
    """

    Parameter('Fas_0'   ,      0) # 3000 corresponds to 50ng/ml Fas(?)
    Parameter('CD95_0'  ,    200) # 200 receptors per cell
    Parameter('FADD_0'  ,  1.0e3) # molecules per cell (arbitrarily assigned)1000
    Parameter('flip_L_0',  1.0e4) # molecules per cell
    Parameter('flip_S_0',  1.0e4) # molecules per cell
    Parameter('proC8_0' ,  2.0e4) # procaspase 8 molecules per cell 20000
    Parameter('C8_0'    ,      0) # active caspase 8 dimers per cell.
    Parameter('Bar_0'   ,  1.0e3) # Bar molecules per cell.

    Initial(Fas(blig=None), Fas_0)       #Fas Ligand
    Initial(CD95(blig=None, bDD=None), CD95_0)     #Fas Receptor (CD95)
    Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0) #FADD
    Initial(flip_L(bDED=None), flip_L_0)   #c-Flip[L]
    Initial(flip_S(bDED=None), flip_S_0)   #c-Flip[S]
    Initial(proC8(bDED=None), proC8_0)    #procaspase 8
    Initial(C8(bC8=None), C8_0)       #caspase 8
    Initial(Bar(bC8=None), Bar_0)     #bifunctional apoptosis regulator
    
    # =========================================
    # CD95 ligation and formation of Secondary Complex rules
    # -----------------------------------------
    #   Fas + CD95 <-> Fas:CD95
    #   Fas:CD95 + FADD <-> Fas:CD95:FADD
    #   Fas:CD95:FADD + proC8 <-> Fas:CD95:FADD:proC8
    
    #   Fas:CD95:FADD:proC8 + proC8 <-> Fas:CD95:FADD:proC8:proC8 -> Fas:CD95:FADD + C8
    #   Fas:CD95:FADD:proC8 + flip_L <-> Fas:CD95:FADD:proC8:flip_L
    #   Fas:CD95:FADD:proC8 + flip_S <-> Fas:CD95:FADD:proC8:flip_S
    
    #   Fas:CD95:FADD:proC8:flip_L <-> Fas:CD95 + FADD:proC8:flip_L
    #   Fas:CD95:FADD:proC8:flip_S <-> Fas:CD95 + FADD:proC8:flip_S
    # ------------------------------------------

    # -------------DISC assembly----------------
    bind(Fas(blig=None), 'blig',  CD95(blig = None, bDD = None), 'blig', [KF, KR])
    bind(CD95(blig = ANY, bDD = None), 'bDD', FADD(bDD = None, bDED1 =None, bDED2 = None), 'bDD', [KF, KR])
    
    bind(FADD(bDD = ANY, bDED1 = None, bDED2 = None),'bDED1', proC8(bDED = None), 'bDED', [KF, KR])
    # For simplicity allow proC8 to bind FADD before any c-Flip do.
    bind(FADD(bDD = ANY, bDED2 = None, bDED1 = ANY), 'bDED2', flip_L(bDED = None), 'bDED', [KF, KR])
    bind(FADD(bDD = ANY, bDED2 = None, bDED1 = ANY), 'bDED2', flip_S(bDED = None), 'bDED', [KF, KR])
    
    
    # procaspase 8 dimerization and activation
    bind(FADD(bDD = ANY, bDED2 = None, bDED1 = ANY), 'bDED2', proC8(bDED = None), 'bDED', [KF, KR])
    DISC_proC8 = CD95(blig=ANY, bDD=ANY) % Fas(blig=ANY) % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % proC8(bDED=ANY)%proC8(bDED=ANY)
    DISC = CD95(blig=ANY, bDD=ANY) % Fas(blig=ANY) % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)
    Rule('C8_activation', FADD(bDED1 = ANY, bDED2 = ANY)%proC8(bDED=ANY)%proC8(bDED = ANY) >> FADD(bDED1 = None, bDED2 = None) + C8(bC8 = None), KC)

    # caspase 8 inhibition by BAR
    bind(Bar(bC8 = None), 'bC8', C8(bC8 = None), 'bC8', [KF, KR])
    
    # release of secondary complex from the DISC
    bind(FADD(bDD = None, bDED2 = ANY, bDED1 = ANY), 'bDD', CD95(blig = ANY, bDD=None), 'bDD', [Parameter('k1', 0),KR])

def TNFR1_to_SecondaryComplex_monomers():
    """ Declares TNFa, TNFR1, TRADD, CompI, RIP1, A20, CYLD, NEMO and NFkB.
    Upon activation, TNFR1 gets endocytosed and post translationally modified After Complex I 
    has released TRADD and RIP1 it possibly gets recycled. This is represented by giving TNFR1 
    two states: norm and spent. TRADD has two states. CompI has two states. RIP1 has three states:
    Ub, PO4 and inactive. RIP1 binds FADD, Complex I, Bid-P and RIP3 (see SecondaryComplex_Bid).
    A20, CYLD and NEMO catalyze transformations of CompI and RIP1, maybe theycan be represented in the rates. 
    """
    Monomer('TNFa', ['blig'])
    Monomer('TNFR1', ['blig', 'bDD', 'state'], {'state':['norm','spent']})
    Monomer('TRADD', ['bDD1', 'bDD2', 'state'], {'state':['active', 'inactive']})
    Monomer('CompI', ['bDD', 'state'], {'state':['unmod', 'mod']}) #Neglecting RIP1:proC8 binding.. for simplicity.
    Monomer('RIP1', ['bDD', 'bRHIM', 'state'], {'state':['unmod', 'ub', 'po4', 'trunc']})
    Monomer('NFkB', ['bf'])

def TNFR1_to_SecondaryComplex():
    """Defines the interactoins from TNFR1 ligation to generation of secondary
    complexes as per ANRM 1.0.
        
    Uses TNFa, TNFR1, TRADD, CompI, RIP1 and NFkB. C8 active dimers and
    their associated parameters to generate rules that describe the ligand/receptor
    binding, FADD recruitment, proC8 and c-flip recruitment, activation of caspase
    and release of the secondary complex. This model assumes that one copy of proC8
    binds FADD before c-flip.
    
    This model converts proC8:proC8 to C8 (active caspase 8 dimer)
    This model also produces Secondary complex, FADD:proC8:c-Flip.
    """
    
    Parameter('TNFa_0'  ,     0) # 3000 corresponds to 50ng/ml TNFa
    Parameter('TNFR1_0' ,   200) # 200 receptors per cell
    Parameter('TRADD_0' ,  1000) # molecules per cell (arbitrarily assigned)1000
    Parameter('CompI_0' ,     0) # complexes per cell
    Parameter('RIP1_0'  , 20000) # molecules per cell 20000
    Parameter('NFkB_0'  ,     0) # molecules per cell
    
    Initial(TNFa(blig=None), TNFa_0)                                 # TNFa Ligand
    Initial(TNFR1(blig=None, bDD=None, state='norm'), TNFR1_0)       # TNFR1
    Initial(TRADD(bDD1=None, bDD2=None, state='inactive'), TRADD_0)  # TRADD
    Initial(CompI(bDD=None, state='unmod'), CompI_0)                 # Complex I
    Initial(RIP1(bDD=None, bRHIM = None, state = 'unmod'), RIP1_0)   # RIP1
    Initial(NFkB(bf=None), NFkB_0)

    # =========================================
    # TNFR1 ligation, formation of Complex I and release of RIP1 and TRADD rules
    # -----------------------------------------
    #   TNFa+ TNFR1 <-> TNFa:TNFR1
    #   TNFa:TNFR1 + TRADD <-> TNFa:TNFR1:TRADD >> CompI
    #   CompI + RIP1 <-> CompI:RIP1 >> [active]CompI:RIP1-Ub
    
    #   [active]CompI:RIP1-Ub >> NFkB # This reaction will consume the receptor.
    #   [active]CompI:RIP1-Ub >> [active]CompI:RIP1
    #   [active]CompI:RIP1 >> [active]CompI # A20 mediated degradation of RIP1
    #   [active]CompI:RIP1 >> [active]CompI + RIP1
    #   [active]CompI >> [active]TRADD + TNFa:[spent]TNFR1
    
    #   TNFa:[spent]TNFR1 >> [norm]TNFR1 #receptor recycle typically distroys the ligand.
    # ------------------------------------------
    
    # -------------Complex I assembly----------------
    bind(TNFa(blig=None), 'blig', TNFR1(blig=None, bDD=None, state='norm'), 'blig', [KF, KR])
    bind(TNFR1(blig = ANY, bDD = None, state =  'norm'), 'bDD', TRADD(bDD1=None, bDD2=None, state='inactive'), 'bDD1', [KF, KR])
    preCompI = TNFa(blig=ANY)%TNFR1(blig=ANY, bDD=ANY, state = 'norm')%TRADD(bDD1 = ANY, bDD2=None, state = 'inactive')
    Rule('CompI_formation', preCompI >> CompI(bDD=None, state = 'unmod'), KC)
    
    # --------------Complex I - RIP1 Modification-----------------
    bind(CompI(bDD=None, state = 'unmod'), 'bDD', RIP1(bDD=None, bRHIM = None, state='unmod'), 'bDD',[KF, KR])
    bind(CompI(bDD=None, state = 'unmod'), 'bDD', RIP1(bDD=None, bRHIM = None, state='ub'), 'bDD',[KF, KR]) 
    
    Rule('CompI_Ub', CompI(bDD=ANY, state = 'unmod')%RIP1(bDD=ANY,bRHIM=None, state = 'unmod')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub'), KC)
    Rule('CompI_Ub2', CompI(bDD=ANY, state = 'unmod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub'), KC)
    Rule('CompI_deUb', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='ub')>>CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,state='unmod'),KC)
    Rule('RIP1_deg', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='unmod') >> CompI(bDD=None, state='mod'),KC2)
    Rule('RIP1_rel', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='unmod') >> CompI(bDD=None, state='mod') + RIP1(bDD=None, bRHIM = None,  state = 'unmod'), KC2)
    Rule('TNFR1_recycle', CompI(bDD=None, state='mod') >> TRADD(bDD1=None, bDD2 = None, state='active') + TNFR1(blig = None, bDD = None, state =  'norm'), KC)
    Rule('NFkB_expression', CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub') + NFkB(bf=None), KE)
    # --------------RIP1 Ubiquitination---------------------------
    Rule('RIP1_Ub', RIP1(bDD=None, bRHIM = None, state='unmod')>> RIP1(bDD=None, bRHIM = None, state='ub'), KC2)
    Rule('RIP1_deUb', RIP1(bDD=None, bRHIM = None, state='ub')>> RIP1(bDD=None, bRHIM = None, state='unmod'), KR)
    
def SecondaryComplex_to_Bid_monomers():
    Monomer('Bid', ['bf', 'state'], {'state':['unmod', 'po4', 'trunc']})
    Monomer('BidK', ['bf']) #unknown Bid-kinase
    Monomer('RIP3', ['bRHIM', 'state'], {'state':['unmod', 'po4', 'trunc']})


def SecondaryComplex_to_Bid():
    """Defines the interactoins from TRADD RIP1 and FADD to generation of secondary
        complexes and truncation of Bid as per ANRM 1.0.
        
        Uses FADD, proC8, C8, flip_L, flip_S, TRADD, RIP1, RIP3 and Bid, and their
        associated parameters to generate rules that describe the FADD recruitment, 
        proC8 and c-flip recruitment to the secondary complex, activation of caspase
        truncation of RIP1 and RIP3 and phosphorylation of Bid. This model assumes 
        that one copy of proC8 binds FADD before c-flip. Among other things...
        
        This model converts proC8:proC8 to C8 (active caspase 8 dimer)
        This model also produces Secondary complex, FADD:proC8:c-Flip.
        """
    Parameter('RIP3_0'  , 2.0e4) # molecules per cell
    Parameter('Bid_0'   , 2.0e4) # molecules per cell
    Parameter('BidK_0'  , 5.0e3) # molecules per cell
    
    Initial(RIP3(bRHIM = None, state = 'unmod'), RIP3_0)   # RIP3
    Initial(Bid(bf = None, state = 'unmod'), Bid_0)        # Bid
    Initial(BidK(bf = None), BidK_0)
    

    # -------------Assembling Complex II-----------------
    bind(FADD(bDD = None, bDED1 = None, bDED2 = None), 'bDD', TRADD(bDD1=None, state = 'active'), 'bDD1', [KF, KR])
    bind(FADD(bDD = None, bDED1 = None, bDED2 = None), 'bDD', RIP1(bDD=None, bRHIM = None, state = 'unmod'), 'bDD', [Ka_RIP1_FADD, Kd_RIP1_FADD])
    bind(TRADD(bDD2 = None, state = 'active'),'bDD2', RIP1(bDD = None, bRHIM = None, state = 'unmod'), 'bDD', [KF, KR])
    # For simplicity, I am neglecting the binary intereaction that occurs between proC8 and RIP1.
    # Binding of proC8 and c-flip to FADD is accomplished in CD95_to_Secondary complex. 

    #--------------RIP1 Truncation reactions-------------
    #---Truncation by C8---------------------------------
    CIIA = TRADD(bDD2 = None, bDD1 = ANY, state = 'active') %  FADD(bDD=ANY, bDED1=None, bDED2=None)
    RIP_CIIA_proC8 = RIP1(bDD=ANY, bRHIM = None, state = 'unmod')% TRADD(bDD2 = None, bDD1 = ANY, state = 'active') % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED=ANY)%proC8(bDED=ANY)
    RIP_CIIB_proC8 = RIP1(bDD=ANY, bRHIM = None, state = 'unmod')% FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED=ANY)%proC8(bDED=ANY)
    Rule('RIP1_truncation_CIIA', RIP_CIIA_proC8 >> CIIA + C8(bC8=None) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), KC)
    Rule('RIP1_truncation_CIIB', RIP_CIIB_proC8 >> FADD(bDD=None, bDED1=None, bDED2=None)+ C8(bC8=None) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), KC)
    
    catalyze_state(C8(bC8=None), 'bC8', RIP1(bDD=None), 'bRHIM', 'state', 'unmod', 'trunc', [KF, KR, KC])

    #---Truncation by proC8:cFlip_L---------------------
    Riptosome_FADD = RIP1(bDD=1, bRHIM = None, state = 'unmod')%FADD(bDD=1, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY)
    Rule('RIP1_truncation2', Riptosome_FADD >> FADD(bDD=None, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), KC)

    Riptosome_TRADD = RIP1(bDD=1, bRHIM = None, state = 'unmod')%TRADD(bDD1=ANY, bDD2=1)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY)
    Rule('RIP1_truncation1', Riptosome_TRADD >> FADD(bDD=None, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), KC)
    
    # -------------RIP3 Binding Interactions
    Ripto1_Flip_S = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=None, state='unmod') % TRADD(bDD1=ANY, bDD2=ANY, state='active') % flip_S(bDED=ANY) % proC8(bDED=ANY)
    Necrosome1 = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=6, state='unmod') % TRADD(bDD1=ANY, bDD2=ANY, state='active') % flip_S(bDED=ANY) % proC8(bDED=ANY) % RIP3(bRHIM= 6, state = 'unmod')
    Rule('RIP3_binding1', Ripto1_Flip_S + RIP3(bRHIM= None, state = 'unmod') <> Necrosome1, KF, KR)

    Ripto2_Flip_S = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=None, state='unmod') % flip_S(bDED=ANY) % proC8(bDED=ANY)
    Necrosome2 = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=5, state='unmod') % flip_S(bDED=ANY) % proC8(bDED=ANY) % RIP3(bRHIM= 5, state = 'unmod')
    Rule('RIP3_binding2', Ripto2_Flip_S + RIP3(bRHIM= None, state = 'unmod') <> Necrosome2, KF, KR)
    
    #RIP3 Truncation
    catalyze_state(C8(bC8=None), 'bC8', RIP3(), 'bRHIM', 'state', 'unmod', 'trunc', [KF, KR, KC])

    # Bid Phosphorylation and Truncation
    catalyze_state(BidK(), 'bf', Bid(), 'bf', 'state', 'unmod', 'po4', [KF, KR, KC])
    catalyze_state(C8(bC8=None), 'bC8', Bid(), 'bf', 'state', 'unmod', 'trunc', [KF, KR, KC])

    # Bid-PO4 competing with RIP1 for binding to Complex II
    bind(TRADD(bDD2 = None, state = 'active'),'bDD2', Bid(bf = None, state = 'po4'), 'bf', [KF, KR])
    # Bid-PO4 sequestering RIP1
    bind(RIP1(bDD = None, bRHIM = None, state = 'unmod'), 'bRHIM', Bid(bf = None, state = 'po4'), 'bf', [KF, KR])

def apaf1_to_parp_monomers():
    """ Declares CytochromeC, Smac, Apaf-1, the Apoptosome, Caspases 3, 6, 9,
    XIAP and PARP.

    The package variable 'bf' specifies the name of the site to be used
    for all binding reactions.

    The 'state' site denotes various localization and/or activity states of a
    Monomer, with 'C' denoting cytoplasmic localization and 'M' mitochondrial
    localization.
    """

    # Cytochrome C
    Monomer('Apaf', ['bf', 'state'], {'state':['I', 'A']}) # Apaf-1
    Monomer('Apop', ['bf']) # Apoptosome (activated Apaf-1 + caspase 9)
    # Csp 3, states: pro, active, ubiquitinated
    Monomer('C3', ['bf', 'state'], {'state':['pro', 'A', 'ub']})
    # Caspase 6, states: pro-, Active
    Monomer('C6', ['bf', 'state'], {'state':['pro', 'A']})
    Monomer('C9', ['bf']) # Caspase 9
    # PARP, states: Uncleaved, Cleaved
    Monomer('PARP', ['bf', 'state'], {'state':['U', 'C']})
    Monomer('XIAP', ['bf']) # X-linked Inhibitor of Apoptosis Protein
    
def pore_to_parp():
    """Defines what happens after the pore is activated and Cytochrome C and
    Smac are released.

    Uses CytoC, Smac, Apaf, Apop, C3, C6, C8, C9, PARP, XIAP monomers and their
    associated parameters to generate the rules that describe apoptosome
    formation, XIAP inhibition, activation of caspases (including
    caspase-6-mediated feedback), and cleavage of effector caspase substrates
    as specified in EARM 1.0.

    Declares initial conditions for CytoC, Smac, Apaf-1, Apoptosome, caspases
    3, 6, and 9, XIAP, and PARP.
    """

    # Declare initial conditions:

    Parameter('Apaf_0'  , 1.0e5) # Apaf-1
    Parameter('C3_0'    , 1.0e4) # procaspase-3 (pro-C3)
    Parameter('C6_0'    , 1.0e4) # procaspase-6 (pro-C6)
    Parameter('C9_0'    , 1.0e5) # procaspase-9 (pro-C9)
    Parameter('XIAP_0'  , 1.0e5) # X-linked inhibitor of apoptosis protein
    Parameter('PARP_0'  , 1.0e6) # C3* substrate

    alias_model_components()

    Initial(Apaf(bf=None, state='I'), Apaf_0)
    Initial(C3(bf=None, state='pro'), C3_0)
    Initial(C6(bf=None, state='pro'), C6_0)
    Initial(C9(bf=None), C9_0)
    Initial(PARP(bf=None, state='U'), PARP_0)
    Initial(XIAP(bf=None), XIAP_0)

    # CytoC and Smac activation after release
    # --------------------------------------

    equilibrate(Smac(bf=None, state='C'), Smac(bf=None, state='A'), [KF, KF])
    equilibrate(CytoC(bf=None, state='C'), CytoC(bf=None, state='A'),[KF, KF])

    # Apoptosome formation
    # --------------------
    #   Apaf + cCytoC <-->  Apaf:cCytoC --> aApaf + cCytoC
    #   aApaf + pC9 <-->  Apop
    #   Apop + pC3 <-->  Apop:pC3 --> Apop + C3

    catalyze_state(CytoC(), 'bf', Apaf(), 'bf', 'state', 'I', 'A', [Kf_Apaf_acti, KR, KC])
    bind(Apaf(bf=None, state='A'),'bf', C9(bf=None), 'bf', [Kf_Apop_asse, KR])
    Rule('Apoptosome', Apaf(bf=ANY, state = 'A')%C9(bf=ANY)>>Apop(bf=None), KC)
    #one_step_conv.. I could not find this macro in the tutorial. I think it removed.
    #one_step_conv(Apaf(state='A'), C9(), Apop(bf=None), [Kf_Apop_asse, KR])
    catalyze_state(Apop(), 'bf', C3(), 'bf', 'state', 'pro', 'A', [Kf_C3_activa, KR, KC]) 

    # Apoptosome-related inhibitors
    # -----------------------------
    #   Apop + XIAP <-->  Apop:XIAP  
    #   cSmac + XIAP <-->  cSmac:XIAP  

    bind(Apop(bf=None), 'bf', XIAP(bf=None), 'bf', [Kf_Apop_inhi, KR]) 
    bind(Smac(bf=None, state='A'), 'bf', XIAP(bf=None), 'bf', [Kf_Smac_inhi, KR]) 

    # Caspase reactions
    # -----------------
    # Includes effectors, inhibitors, and feedback initiators:
    #
    #   pC3 + C8 <--> pC3:C8 --> C3 + C8 CSPS
    #   pC6 + C3 <--> pC6:C3 --> C6 + C3 CSPS
    #   XIAP + C3 <--> XIAP:C3 --> XIAP + C3_U CSPS
    #   PARP + C3 <--> PARP:C3 --> CPARP + C3 CSPS
    #   pC8 + C6 <--> pC8:C6 --> C8 + C6 CSPS
    catalyze_state(C8(), 'bC8', C3(),'bf', 'state', 'pro','A', [Kf_C3_activ2, KR, KC])
    catalyze_state(XIAP(), 'bf', C3(), 'bf', 'state', 'A', 'ub', [Kf_C3_ubiqui, KR, Kc_C3_ubiqui])
    catalyze_state(C3(state='A'), 'bf', PARP(),'bf', 'state', 'U', 'C', [KF, Kr_PARP_clea, KC])
    catalyze_state(C3(state='A'), 'bf', C6(), 'bf', 'state', 'pro', 'A', [KF, KR, KC])
    #catalyze(C6(state='A'), C8(state='pro'), C8(state='A'), [Kf_C8_activ2, KR, KC])

def momp_monomers():
    """Declare the monomers used in the Albeck MOMP modules."""

    # == Activators
    # Bid, states: Untruncated, Truncated, truncated and Mitochondrial
    # Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']})
    # == Effectors
    # Bax, states: Cytoplasmic, Mitochondrial, Active
    # sites 's1' and 's2' are used for pore formation
    Monomer('Bax', ['bf', 's1', 's2', 'state'], {'state':['C', 'M', 'A']})
    # == Anti-Apoptotics
    Monomer('Bcl2', ['bf'])

    # == Cytochrome C and Smac
    Monomer('CytoC', ['bf', 'state'], {'state':['M', 'C', 'A']})
    Monomer('Smac', ['bf', 'state'], {'state':['M', 'C', 'A']})
    
# MOMP module implementations
# ===========================

# Motifs
# ------

# Because several of the models in [Albeck2008]_ overlap, some mechanistic
# aspects have been refactored into the following "motifs", implemented as
# functions:

#def Bax_tetramerizes(bax_active_state='A', rate_scaling_factor=1):
def Bax_tetramerizes():
    """Creates rules for the rxns Bax + Bax <> Bax2, and Bax2 + Bax2 <> Bax4.

    Parameters
    ----------
    bax_active_state : string: 'A' or 'M'
        The state value that should be assigned to the site "state" for
        dimerization and tetramerization to occur.
    rate_scaling_factor : number
        A scaling factor applied to the forward rate constants for dimerization
        and tetramerization. 
    """
    rate_scaling_factor = 1
    active_unbound = {'state': bax_active_state, 'bf': None}
    active_bax_monomer = Bax(s1=None, s2=None, bf = None, state = 'A')
    bax2 =(Bax(s1=1, s2=None, bf = None, state = 'A') %
           Bax(s1=None, s2=2, bf = None, state = 'A'))
    bax4 =(Bax(s1=1, s2=4, bf = None, state = 'A') %
           Bax(s1=2, s2=1, bf = None, state = 'A') %
           Bax(s1=3, s2=2, bf = None, state = 'A') %
           Bax(s1=4, s2=3, bf = None, state = 'A'))
    
    #active_bax_monomer = Bax(s1=None, s2=None, **active_unbound)
    #bax2 =(Bax(s1=1, s2=None, **active_unbound) %
    #       Bax(s1=None, s2=1, **active_unbound))
    #bax4 =(Bax(s1=1, s2=4, **active_unbound) %
    #       Bax(s1=2, s2=1, **active_unbound) %
    #       Bax(s1=3, s2=2, **active_unbound) %
    #       Bax(s1=4, s2=3, **active_unbound))
    KF_scaled = 1e-6*rate_scaling_factor
    Rule('Bax_dimerization', active_bax_monomer + active_bax_monomer <> bax2,
         Parameter('Bax_dimerization_kf', KF_scaled),
         Parameter('Bax_dimerization_kr', KR))
    # Notes on the parameter values used below:
    #  - The factor 2 is applied to the forward tetramerization rate because
    #    BNG (correctly) divides the provided forward rate constant by 1/2 to
    #    account for the fact that Bax2 + Bax2 is a homodimerization reaction,
    #    and hence the effective rate is half that of an analogous
    #    heterodimerization reaction. However, Albeck et al. used the same
    #    default rate constant of 1e-6 for this reaction as well, therefore it
    #    must be multiplied by 2 in order to match the original model
    #  - BNG apparently applies a scaling factor of 2 to the reverse reaction
    #    rate, for reasons we do not entirely understand. The factor of 0.5 is
    #    applied here to make the rate match the original Albeck ODEs.
    KF_Bax_tetramerization = 2*1.0e-6*rate_scaling_factor
    KR_Bax_tetramerization = 0.5*1.0e-3
    Rule('Bax_tetramerization', bax2 + bax2 <> bax4,
         Parameter('Bax_tetramerization_kf', KF_Bax_tetramerization),
         Parameter('Bax_tetramerization_kr', KR_Bax_tetramerization))

def Bcl2_binds_Bax1_Bax2_and_Bax4(bax_active_state='A', rate_scaling_factor=1):
    """Creates rules for binding of Bcl2 to Bax monomers and oligomers.

    Parameters
    ----------
    bax_active_state : string: 'A' or 'M'
        The state value that should be assigned to the site "state" for
        the Bax subunits in the pore.
    rate_scaling_factor : number
        A scaling factor applied to the forward rate constants for binding
        between Bax (monomers, oligomers) and Bcl2.
    """
    KF_scaled = 1.0e-6*rate_scaling_factor
    Parameter('KF_Bcl2binding', KF_scaled)
    bind(Bax(state=bax_active_state, s1=None, s2=None), 'bf', Bcl2, 'bf',
         [KF_Bcl2binding, KR])
    pore_bind(Bax(state=bax_active_state), 's1', 's2', 'bf', 2, Bcl2, 'bf',
         [KF_Bcl2binding, KR])
    pore_bind(Bax(state=bax_active_state), 's1', 's2', 'bf', 4, Bcl2, 'bf',
         [KF_Bcl2binding, KR])

def albeck_11c(do_pore_transport=True):
    """Model incorporating Bax oligomerization.

    Features:
        - Bid activates Bax
        - Active Bax dimerizes; Bax dimers dimerize to form tetramers
        - Bcl2 binds/inhibits Bax monomers, dimers, and tetramers
        - Bax tetramers bind to and transport Smac to the cytosol
    """

    #alias_model_components()
    #Initial(Bid(state='U', bf=None), Parameter('Bid_0', 4e4))
    Bax_Inactive = 'C'
    inactive_monomer = {'state': Bax_Inactive, 'bf': None}
    Initial(Bax(s1 = None, s2 = None, **inactive_monomer), Parameter('Bax_0', 1e5))
    Initial(Bcl2(bf=None), Parameter('Bcl2_0', 2e4))

    # tBid activates Bax
    catalyze_state(Bid(state='trunc'),'bf', Bax(),'bf', 'state', 'C', 'A',
             [Kf_Bax_activ, KR, KC])

    # Bax dimerizes/tetramerizes
    Bax_tetramerizes(bax_active_state='A')

    # Bcl2 inhibits Bax, Bax2, and Bax4
    Bcl2_binds_Bax1_Bax2_and_Bax4(bax_active_state='A')
    
    if do_pore_transport:
        Initial(Smac(state='M', bf=None), Parameter('Smac_0', 1e6))
        Initial(CytoC(state='M', bf=None), Parameter('CytoC_0', 1e6))
        
        # NOTE change in KF here from previous model!!!!
        pore_transport(Bax(state='A'), 's1','s2','bf', 4, 4, Smac(state='M'), 'bf', Smac(state='C'),
            [[Parameter('KF_Smac_transport', 2.0e-6), KR, Parameter('KC_Smac_transport', 10)]])
        pore_transport(Bax(state='A'), 's1','s2','bf', 4, 4, CytoC(state='M'), 'bf', CytoC(state='C'),
            [[KF, KR, KC_Smac_transport]])
        
CD95_to_SecondaryComplex_monomers()
CD95_to_SecondaryComplex()
TNFR1_to_SecondaryComplex_monomers()
TNFR1_to_SecondaryComplex()

SecondaryComplex_to_Bid_monomers()
SecondaryComplex_to_Bid()

momp_monomers()
apaf1_to_parp_monomers()
pore_to_parp()
#albeck_11c()
Bax_tetramerizes()
Bcl2_binds_Bax1_Bax2_and_Bax4(bax_active_state='A')

Observable('Obs_TNFa', TNFa(blig =  None))
Observable('Obs_Fas', Fas(blig = None))
Observable('Obs_TNFR1', TNFR1(blig = None))
Observable('Obs_CD95', CD95(blig = None))
Observable('CD95_Fas', CD95(blig = ANY))
Observable('DISC', CD95(blig = ANY, bDD=ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2 = ANY))
Observable('TNFR1_TNF', TNFR1(blig=ANY,bDD = ANY))
Observable('ComplexI', CompI())
Observable('Obs_RIP1', RIP1(state = 'unmod'))
Observable('Obs_TRADD', TRADD(state = 'inactive'))
Observable('Obs_TRADDa', TRADD(state = 'active'))
Observable('SecondaryComplex', FADD(bDD=None, bDED1 = ANY, bDED2 = ANY))
Observable('Complex_IIA', TRADD(bDD1=ANY, bDD2=None)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
Observable('Riptosome1', RIP1(bDD = ANY, bRHIM = None)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
Observable('Riptosome2', RIP1(bDD = ANY, bRHIM = None)%TRADD(bDD1=ANY, bDD2=ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
Observable('Obs_RIP1_Ub', RIP1(state = 'ub'))
Observable('CompI_RIP1', CompI(bDD=ANY))
Observable('CompI_mod', CompI(state='mod'))
Observable('RIP1_Bid', RIP1()%Bid())
Observable('Bid_Riptosome1', Bid(bf= ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
Observable('Bid_Riptosome2', Bid(bf= ANY)%TRADD(bDD1=ANY, bDD2=ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
Observable('Bid_Trunc', Bid(state='trunc'))
Observable('Bid_PO4', Bid(state='po4'))
Observable('RIP1_Trunc', RIP1(state='trunc'))
Observable('RIP3_Trunc', RIP3(state='trunc'))
Observable('Obs_proC8', proC8())
Observable('Obs_C8', C8())
Observable('FADD_proC8_proC8', FADD(bDED1 = ANY, bDED2 = ANY)%proC8(bDED=ANY)%proC8(bDED = ANY))
