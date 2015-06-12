#!/bin/sh

#################################################################################
#                                                                               #
#   get_Newbler_gene_to_trans_map.sh                                            #
#                                                                               #
#   Takes as input a 454Isotigs.txt resulting from 454's gsAssembler and makes  #
#   the isoform - gene map in the same way as it is done with Trinity and       #
#   Trinotate.                                                                  #
#                                                                               #
#################################################################################

fileIn=$1

paste \
    <(  cat $fileIn                 \
        | grep ^">"                 \
        | awk ' { print $2 } '      \
        | sed ' s/gene=//    '    ) \
    <(  cat $fileIn                 \
        | grep  ^">"                \
        | awk   ' { print $1 } '    \
        | sed   ' s/>// '         )

