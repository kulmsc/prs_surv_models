author=$1

#Rscript get_y_data.R $author
#python zero_enet_cox.py
#Rscript analysis.R $author
#mv static_pred.txt.gz res/${author}.enet_pred.gz
#mv coef.txt.gz res/${author}.coef.gz
#mv train_pred.txt.gz res/${author}.train_pred.gz
#mv best_alpha.txt res/${author}.best_alpha.txt

#go backwards
Rscript get_y_data.R $author
cp res/${author}.enet_pred.gz static_pred.txt.gz
cp res/${author}.coef.gz coef.txt.gz
cp res/${author}.train_pred.gz train_pred.txt.gz
#Rscript analysis.R $author
Rscript gof.R $author
#Rscript  attach_pheno.R $author
#rm coef.txt.gz static_pred.txt.gz train_pred.txt.gz
