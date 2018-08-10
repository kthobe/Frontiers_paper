
# Script for the generation and classification of a model pool
# This code was used to generate the results for the Paper "Evaluating Uncertainty in Signaling Networks using Logical Modeling" by Thobe et al., 2018

# TomClass was developed by Hannes Klarner: https://github.com/hklarner/TomClass
# for more background information see the Manual of TomClass: https://github.com/hklarner/TomClass/blob/master/Documentation/Manual.pdf
# theoretical background is given in Klarner, H. (2014). Contributions to the Analysis of Qualitative Models of Regulatory Networks. Ph.D. thesis, Freie Universitaet Berlin
# as well as Thobe, K. (2017). Logical modeling of uncertainty in signaling pathways of cancer systems. Ph.D. thesis,Freie Universitaet Berlin

# This script was edited by Kirsten Thobe, 01.08.2018

import sys
import itertools as IT
import random as RND

import Instantiation
import Annotation
import Analysis
from Engine import Database
from Engine import ModelChecking as MCheck
import Engine
import time


UPDATE = 'asynchronous'                 # in current TomClass version only asynchronous possible

DBNAME = 'paper.sqlite'                 # enter SQL file for storing data base
FNAMECSV = 'full.csv'                   # enter CVS file for storing classification table
                                        # full.cvs for no RESTRICTION
                                        # 1257.cvs for cell line MZ1257RC
                                        # 1851.cvs for cell line MZ1851RC

create                  = 0             # select 1 for creating the database, only necessary once

annotate_compatible     = 0             # select for model checking of data on pool
annotate_crosstalk      = 0             # select for annotating presence/absence of optional edges to database, prerequisite for classification


analyse_classes         = 1             # select for classification analysis, enter classes under CLASSES
export_csv              = 1             # select for creating CVS file

# For creating classifications for a subset of model, select classes (e.g. presence of edges, data sets) as RESTRICTION
# for cell line MZ1257RC enter 'WBDMSO and WB1257Sora and Bp1257Sora and Bp1257DMSO and Bp1257Sora2 and Bp1257DMSO2'
# for cell line MZ1851RC enter 'WBDMSO and WB1851Sora and Bp1851Sora and Bp1851DMSO and Bp1851Sora2 and Bp1851DMSO2'
RESTRICTION             = ''

# List all annotated properties (data sets, edges) the model pool is going to be classified for
# Example:              = ['A', 'B', 'C']
CLASSES                 = ['Crosstalk','Ct_Sora_Raf','Ct_Sora_EGFR','Ct_Sora_IGFR','Ct_EGFR_PI3K','Ct_Erk_mTor','Ct_mTor_IGFR','Ct_PI3K_Raf',
'WBDMSO','WB1257Sora','WB1851Sora','Bp1851Sora','Bp1851DMSO','Bp1851Sora2','Bp1851DMSO2','Bp1257Sora', 'Bp1257DMSO','Bp1257Sora2', 'Bp1257DMSO2']

# program for creating a CTL formula from dataset--------------------------------------------------------------
# creates a CTL formula from a matrix, do not change
def BooleanTSformula( Names, Measurements, Fixpoint ):
    Monotone=False

    for i,m in enumerate(Measurements):
        if not len(Names)==len(m):
            print 'Measurement %i is of length %i but %i names are given.'%(i,len(m),len(Names))
            raise Exception
    if len(Measurements)<2:
        print 'Need at least 2 measurements.'
        raise Exception

    s = '?'
    for i,m in enumerate(Measurements):
        m_str = '&'.join(['%s=%i'%item for item in zip(Names,m)])
        assert('?' not in m_str)
        assert(s.count('?')==1)

        if i==0:
            init = m_str
        else:
            if Monotone:
                cond = '&'.join(  ['(%s=%i->X %s=%i)'%(Names[j],v,Names[j],v) for j,v in enumerate(m)]  )
                rep  = '((%s) U (%s?))'%(cond,m_str)
                if i>1: rep='&'+rep
                s=s.replace('?',rep)
            else:
                rep = 'EF(%s?)'%m_str
                if i>1: rep='&'+rep
                s=s.replace('?', rep)


    assert(s.count('?')==1)
    if Fixpoint:
    	s=s.replace('?','&Delta=0')
    else:
    	s=s.replace('?','')

    return init, s

    #-----------------------------------------------------------------------------------------------------------

def run():
    if create:
        # Model defintion, first components and number of states (1 for Boolean)
        parameters = {'Db_name'         : DBNAME,
                      'Components'      : [('EGF',1),
                                           ('EGFR',1),
                                           ('IGF',1),
                                           ('IGFR',1),
                                           ('Sora',1),
                                           ('Raf',1),
                                           ('Erk',1),
                                           ('PI3K',1),
                                           ('Akt',1),
                                           ('mTor',1)],
        # definition of interactions with threshold
                      'Interactions'    : [('EGF','EGF',(1,)),
                                           ('EGF','EGFR',(1,)),
                                           ('IGF','IGF',(1,)),
                                           ('IGF','IGFR',(1,)),
                                           ('Sora','Sora',(1,)),
                                           ('EGFR','Raf',(1,)),
                                           ('Raf','Erk',(1,)),
                                           ('Erk','EGFR',(1,)),
                                           ('IGFR','PI3K',(1,)),
                                           ('PI3K','Akt',(1,)),
                                           ('Akt','mTor',(1,)),

                                           # Crosstalk
                                           ('mTor','IGFR',(1,)),
                                           ('PI3K','Raf',(1,)),
                                           ('Sora','Raf',(1,)),
                                           ('Sora','EGFR',(1,)),
                                           ('Sora','IGFR',(1,)),
                                           ('EGFR','PI3K',(1,)),
                                           ('Erk','mTor',(1,))
                                           ]}

        # Defintion of edge labels and logical functions:
        # Boolean(A=1,B) states that B has only A as regulator and becomes active for A=1
        # ActivatingOnly and InhibitingOnly defines essential edges, but does not define the logical function
        # NotActivating and NotInhibiting creates optional edge labels
        clauses = ['Boolean(EGF=1, EGF)',
                   'Boolean(IGF=1, IGF)',
                   'Boolean(Sora=1, Sora)',
                   'Boolean(Raf=1, Erk)',
                   'Boolean(PI3K=1, Akt)',
                   'ActivatingOnly(IGF,IGFR,1)',
                   'ActivatingOnly(EGF,EGFR,1)',
                   'ActivatingOnly(EGFR,Raf,1)',
                   'ActivatingOnly(IGFR,PI3K,1)',
                   'InhibitingOnly(Erk,EGFR,1)',
                   'ActivatingOnly(Akt,mTor,1)',
                   # Crosstalk
                   'NotActivating(mTor,IGFR,1)',
                   'NotActivating(Sora,EGFR,1)',
                   'NotActivating(Sora,IGFR,1)',
                   'NotActivating(Sora,Raf,1)',
                   'NotInhibiting(EGFR,PI3K,1)',
                   'NotInhibiting(Erk,mTor,1)',
                   'NotInhibiting(PI3K,Raf,1)'
 ]

        parameters['Constraint'] = ' and '.join( clauses )

        Instantiation.CPEnumeration.run( parameters )

    if not create:
        db = Database.Interface(DBNAME)
        model = db.get_model()
        model.info()
        db.close()

    if annotate_compatible:
    # this function performs model checking for CTL formulas on the model pool
    # for detailed description, see Manual

    # annotate steady-state data -------------------------------------------

        parameters = {  'Db_name'           : DBNAME,
                        'Restriction'       : RESTRICTION,
                        'Property_name'      : 'WBDMSO',
                        'Description'       : 'mTor measurement for DMSO.',
                        'Formula'           : 'EF( Delta=0& mTor=1)',
                        'Initial_states'     : 'Sora=0',
                        'Verification_type'  : 'forsome',
                        'Fix'               : {}}
        Annotation.CTL.run( parameters )


        parameters = {  'Db_name'           : DBNAME,
                        'Restriction'       : RESTRICTION,
                        'Property_name'      : 'WB1851Sora',
                        'Description'       : 'mTor measurement for 1851 Sora.',
                        'Formula'           : 'EF( Delta=0&mTor=1)',
                        'Initial_states'     : 'Sora=1',
                        'Verification_type'  : 'forsome',
                        'Fix'               : {}}
        Annotation.CTL.run( parameters )

        parameters = {  'Db_name'           : DBNAME,
                        'Restriction'       : RESTRICTION,
                        'Property_name'      : 'WB1257Sora',
                        'Description'       : 'mTor measurement for 1257 Sora.',
                        'Formula'           : 'EF( Delta=0& mTor=0)',
                        'Initial_states'     : 'Sora=1',
                        'Verification_type'  : 'forsome',
                        'Fix'               : {}}
        Annotation.CTL.run( parameters )

        parameters = {  'Db_name'           : DBNAME,
                        'Restriction'       : RESTRICTION,
                        'Property_name'      : 'Bp1257DMSO',
                        'Description'       : 'Bioplex measurement for 1257 DMSO.',
                        'Formula'           : 'EF( Delta=0&mTor=1&Akt=1&EGFR=1&Erk=1)',
                        'Initial_states'     : 'Sora=0',
                        'Verification_type'  : 'forsome',
                        'Fix'               : {}}
        Annotation.CTL.run( parameters )

    # annotate time-series data -------------------------------------------
    # 1851 Experiment 1
        names        = [ 'Erk', 'EGFR', 'mTor', 'Akt', 'IGFR']
        measurements = [[1,      1,      1,       1,     1],
                        [1,      1,      1,       1,     0],
                        [1,      1,      1,       0,     0],
                        [0,      0,      0,       0,     0],
                        ]
        init, formula = BooleanTSformula( names, measurements, Fixpoint=False)
        print init, formula

        parameters = {  'Db_name'           : DBNAME,
                        'Restriction'       : RESTRICTION,
                        'Property_name'      : 'Bp1851DMSO',
                        'Description'       : 'Bioplex measurement for 1851 DMSO.',
                        'Formula'           : formula,
                        'Initial_states'     : init,
                        'Verification_type'  : 'forsome',
                        'Fix'               : {'Sora':0}}
        Annotation.CTL.run( parameters )

        names        = [ 'Erk', 'EGFR', 'mTor', 'Akt', 'IGFR']
        measurements = [[0,      0,      1,       0,     1],
                        [1,      1,      1,       1,     1],
                        [0,      0,      0,       0,     1],
                        ]
        init, formula = BooleanTSformula( names, measurements, Fixpoint=False)
        print init, formula

        parameters = {  'Db_name'           : DBNAME,
                        'Restriction'       : RESTRICTION,
                        'Property_name'      : 'Bp1851Sora',
                        'Description'       : 'Bioplex measurement for 1851 Sora.',
                        'Formula'           : formula,
                        'Initial_states'     : init,
                        'Verification_type'  : 'forsome',
                        'Fix'               : {'Sora':1}}
        Annotation.CTL.run( parameters )

    # 1851 Experiment 2
        names        = [ 'Erk', 'EGFR', 'mTor', 'Akt', 'IGFR']
        measurements = [[1,      1,      1,       1,     0],
                        [1,      0,      0,       1,     0],
                        [1,      0,      0,       0,     0],
                        [1,      1,      1,       1,     0],
                        [0,      0,      0,       0,     0],
                        [1,      1,      0,       0,     0],
                        [1,      1,      1,       1,     0],
                        ]
        init, formula = BooleanTSformula( names, measurements, Fixpoint=False)
        print init, formula

        parameters = {  'Db_name'           : DBNAME,
                        'Restriction'       : RESTRICTION,
                        'Property_name'      : 'Bp1851DMSO2',
                        'Description'       : 'Bioplex measurement for 1851 DMSO.',
                        'Formula'           : formula,
                        'Initial_states'     : init,
                        'Verification_type'  : 'forsome',
                        'Fix'               : {'Sora':0}}
        Annotation.CTL.run( parameters )

        names        = [ 'Erk', 'EGFR', 'mTor', 'Akt', 'IGFR']
        measurements = [[1,      1,      1,       1,     1],
                        [1,      1,      0,       1,     0],
                        [1,      1,      1,       1,     0],
                        [0,      0,      1,       0,     0],
                        [1,      1,      0,       0,     0],
                        [1,      1,      1,       0,     1],
                        ]
        init, formula = BooleanTSformula( names, measurements, Fixpoint=False)
        print init, formula

        parameters = {  'Db_name'           : DBNAME,
                        'Restriction'       : RESTRICTION,
                        'Property_name'      : 'Bp1851Sora2',
                        'Description'       : 'Bioplex measurement for 1851 Sora.',
                        'Formula'           : formula,
                        'Initial_states'     : init,
                        'Verification_type'  : 'forsome',
                        'Fix'               : {'Sora':1}}
        Annotation.CTL.run( parameters )

    # 1257 Experiment 1
        names        = [ 'Erk', 'EGFR', 'mTor', 'Akt']
        measurements = [[1,      1,      1,       1],
                        [0,      0,      0,       0],
                        [0,      1,      0,       0],
                        [1,      1,      1,       1],
                        [0,      1,      0,       0],
                        [1,      1,      1,       0],
                        ]
        init, formula = BooleanTSformula( names, measurements, Fixpoint=False)
        print init, formula

        parameters = {  'Db_name'           : DBNAME,
                        'Restriction'       : RESTRICTION,
                        'Property_name'      : 'Bp1257Sora',
                        'Description'       : 'Bioplex measurement for 1257 Sora.',
                        'Formula'           : formula,
                        'Initial_states'     : init,
                        'Verification_type'  : 'forsome',
                        'Fix'               : {'Sora':1}}
        Annotation.CTL.run( parameters )

    # 1257 Experiment 2
        names        = [ 'Erk', 'EGFR', 'mTor', 'Akt']
        measurements = [[1,      1,      0,       0],
                        [0,      1,      1,       0],
                        [1,      1,      1,       1],
                        [0,      0,      0,       0],
                        [0,      1,      0,       0],
                        [1,      1,      1,       1],
                        ]
        init, formula = BooleanTSformula( names, measurements, Fixpoint=False)
        print init, formula

        parameters = {  'Db_name'           : DBNAME,
                        'Restriction'       : RESTRICTION,
                        'Property_name'      : 'Bp1257DMSO2',
                        'Description'       : 'Bioplex measurement for 1257 DMSO.',
                        'Formula'           : formula,
                        'Initial_states'     : init,
                        'Verification_type'  : 'forsome',
                        'Fix'               : {'Sora':0}}
        Annotation.CTL.run( parameters )

        names        = [ 'Erk', 'EGFR', 'mTor', 'Akt']
        measurements = [[0,      0,      0,       0],
                        [0,      0,      0,       1],
                        [1,      1,      1,       1],
                        [0,      1,      1,       1],
                        [1,      1,      1,       1]
                        ]
        init, formula = BooleanTSformula( names, measurements, Fixpoint=False)
        print init, formula

        parameters = {  'Db_name'           : DBNAME,
                        'Restriction'       : RESTRICTION,
                        'Property_name'      : 'Bp1257Sora2',
                        'Description'       : 'Bioplex measurement for 1257 Sora.',
                        'Formula'           : formula,
                        'Initial_states'     : init,
                        'Verification_type'  : 'forsome',
                        'Fix'               : {'Sora':1}}
        Annotation.CTL.run( parameters )

# Annotate presence and absence of optional egdes to models
    if annotate_crosstalk:

        db = Database.Interface( DBNAME )
        model = db.get_model()
        db.close()

        constraints = []
        # here, the optional edges need to be entered in the form of [('A','B'), ...]
        for crosstalk in [('Sora','Raf'),('Sora','EGFR'),('Sora','IGFR'),('EGFR','PI3K'),('PI3K','Raf'),('Erk','mTor'),('mTor','IGFR')]:
            con = Engine.Constraints.Parser.PredicateFormula.parseString( 'Observable(%s,%s,1)'%crosstalk, parseAll=True)[0]
            con.initialize( model )
            constraints.append( con )

        def custom_algorithm( Model, Parametrization, Labels, States ):
            return {'Crosstalk': len([1 for con in constraints if con(Parametrization)])}


        parameters = { 'Db_name'        : DBNAME,
                       'Property_name'  : 'Crosstalk',
                       'Property_type'  : 'int',
                       'Description'    : 'Number of observable crosstalk interactions between the MAPK and mTor pathways.',
                       'Restriction'    : ''}
        Annotation.CustomLoop.run( parameters, custom_algorithm )

        # here, the optional edges need to be entered in the form of [('A','B'), ...], must be identical to the one above
        for crosstalk in [('Sora','Raf'),('Sora','EGFR'),('Sora','IGFR'),('EGFR','PI3K'),('PI3K','Raf'),('Erk','mTor'),('mTor','IGFR')]:
            parameters = {'Db_name'       : DBNAME,
                          'Restriction'   : '',
                          'Property_name' : 'Ct_%s_%s'%crosstalk,
                          'Specification' : 'Observable(%s,%s,1)'%crosstalk,
                          'Description'   : 'Determines whether the crosstalk "%s,%s" is observable.'%crosstalk}
            Annotation.Predicate.run( parameters )

    if analyse_classes:
    # this function performs the classification
        parameters = { 'Db_name'        : DBNAME,
                       'Properties'     : CLASSES,
                       'Restriction'    : RESTRICTION}
        if export_csv:
            parameters['FileName'] = FNAMECSV

        Analysis.Classes.run( parameters )


if __name__=='__main__':
    run()
