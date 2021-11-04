author=$1

Rscript get_y_data.R $author
python zero_enet_cox.py
Rscript analysis.R $author
mv static_pred.txt.gz res/${author}.enet_pred.gz
mv coef.txt.gz res/${author}.coef.gz


#go backwards
#author=tsoi
#cp res/${author}.enet_pred.gz static_pred.txt.gz
#cp res/${author}.coef.gz coef.txt.gz
