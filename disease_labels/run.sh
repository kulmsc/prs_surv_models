author=$1

Rscript get_y_data.R $author
#python zero_enet_cox.py
cp res/${author}.enet_pred.gz static_pred.txt.gz
Rscript do_removal_analysis.R $author
#mv static_pred.txt.gz res/${author}.enet_pred.gz
rm static_pred.txt.gz
