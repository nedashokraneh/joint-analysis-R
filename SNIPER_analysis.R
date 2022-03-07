x <- SNIPER_embedding$kmean_compartment - c(0,SNIPER_embedding$kmean_compartment[-length(SNIPER_embedding$kmean_compartment)])
sum(x!=0)

t = t/100000
best = apply(t,MARGIN=1, FUN=max)
total = apply(t,MARGIN=1, FUN=sum)
best/total
