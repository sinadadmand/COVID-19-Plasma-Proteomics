# from the whole proteome length of N(around 20,000) our group found A (111) differentially expressed proteins, while COMBAT group found B (154) proteins. 
# The intersection between the two is m (57).
# There are 2 two methods to test the significancy of this overlap. First is the hypergeometric way to calculate the p-value which is:

N = 20000
a = 111
b = 153
m = 57

sum(dhyper(m:min(a,b), max(a,b), N - max(a,b), min(a,b)))
# or:

phyper(min(a,b), max(a,b), N - max(a,b), min(a,b)) - phyper(m-1, max(a,b), N - max(a,b), min(a,b))



# Fisher's exact test:

fisher.test(matrix(c(N-a-b+m, a-m, b-m, m), nrow=2), alternative="greater")
# Jaccard_index is:

m / (a+b-m)


# Both methods finds the same p-value. The Jaccard index is also reported as a measurement of similarity between two sets. (number of intersections / number of unions).
