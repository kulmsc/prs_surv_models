author=$1

Rscript get_y_data.R $author
python zero_enet_cox.py
Rscript do_removal_analysis.R $author
mv static_pred.txt.gz res/${author}.enet_pred.gz
