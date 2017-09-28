#!/bin/bash

IFS=' '
complete -W "qualSta trim alignRead sam2bam bamSta assemble assemSta getUnalnCtg rmRedundant pTpG geneCov geneExist subSample gFamExist bam2bed fastaSta sim" eupan

complete -W "qualSta mergeQualSta trim alignRead sam2bam bamSta assemble assemSta mergeAssemSta getUnalnCtg mergeUnalnCtg rmRedundant pTpG geneCov mergeGeneCov geneExist subSample gFamExist bam2bed fastaSta sim" eupanLSF

complete -W "qualSta mergeQualSta trim alignRead sam2bam bamSta assemble assemSta mergeAssemSta getUnalnCtg mergeUnalnCtg rmRedundant pTpG geneCov mergeGeneCov geneExist subSample gFamExist bam2bed fastaSta sim" eupanSLURM


