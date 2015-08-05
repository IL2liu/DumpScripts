#gem_old(544) vs gem_new(1.366)
#read of length = 36 Nucleotids

# samples #208
#total reads in new = 4920506
#total reads in old = 11468981
#overlap = 4845632

require(VennDiagram)
venn.diagram(list(B = 1:11468981, A = 4920506:6623349),fill = c("red", "green"),
             alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, 
             fontfamily =3, filename = "trial.png");

require(venneuler)
v <- venneuler(c(NewGem=4920506, OldGem=11468981, "NewGem&OldGem"=4845632))
plot(v)
