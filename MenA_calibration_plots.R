



mena.counts <- function(onesim.array, age.week.min, age.week.max) {

age.y.min <- ceiling((age.week.min - 1)/12)
age.y.max <- ceiling((age.week.max)/12)-1
subpop <- onesim.array[age.week.min:age.week.max,,]
subpop.tot <- apply(subpop, c(2,3), sum)
subpop.grandtot <- apply(subpop, 3, sum)
x.vals <- as.Date(as.numeric(dimnames(onesim.array)[[3]]), origin="1970-01-01")

pdf(paste(outputdir, "/", vacc_program, " ", vacc_subprogram, " Ages ", age.y.min, "-", age.y.max, ".pdf", sep=""))
plot(x.vals, subpop.tot[1,]/subpop.grandtot, ylim=c(0, max(subpop.tot[1,]/subpop.grandtot, subpop.tot[2,]/subpop.grandtot)),  xlab="Year", ylab="Fraction of population", main="No Protection", pch="")
lines(x.vals, subpop.tot[1,]/subpop.grandtot)
lines(x.vals, subpop.tot[2,]/subpop.grandtot, col="red")
legend("topleft", c("Susceptible", "Colonized"), lty=c(1,1), col=c("black", "red"))

plot(x.vals, subpop.tot[3,]/subpop.grandtot, ylim=c(0, max(subpop.tot[3,]/subpop.grandtot, subpop.tot[4,]/subpop.grandtot)),  xlab="Year", ylab="Fraction of population", main="Low Protection", pch="")
lines(x.vals, subpop.tot[3,]/subpop.grandtot)
lines(x.vals, subpop.tot[4,]/subpop.grandtot, col="red")
legend("topleft", c("Susceptible", "Colonized"), lty=c(1,1), col=c("black", "red"))

plot(x.vals, subpop.tot[5,]/subpop.grandtot, ylim=c(0, max(subpop.tot[5,]/subpop.grandtot, subpop.tot[6,]/subpop.grandtot)),  xlab="Year", ylab="Fraction of population", main="High Protection", pch="")
lines(x.vals, subpop.tot[5,]/subpop.grandtot)
lines(x.vals, subpop.tot[6,]/subpop.grandtot, col="red")
legend("topleft", c("Susceptible", "Colonized"), lty=c(1,1), col=c("black", "red"))

plot(x.vals, subpop.tot[8,]/subpop.grandtot, ylim=c(0, max(subpop.tot[8,]/subpop.grandtot, subpop.tot[9,]/subpop.grandtot)),  xlab="Year", ylab="Fraction of population", main="Vaccinated", pch="")
lines(x.vals, subpop.tot[8,]/subpop.grandtot)
lines(x.vals, subpop.tot[9,]/subpop.grandtot, col="red")
legend("topleft", c("Susceptible", "Colonized"), lty=c(1,1), col=c("black", "red"))

plot(x.vals, subpop.tot[7,]/subpop.grandtot, xlab="Year", ylab="Fraction of population", main="Diseased", pch="")
lines(x.vals, subpop.tot[7,]/subpop.grandtot)
lines(x.vals, subpop.tot[10,]/subpop.grandtot, col="red")
legend("topleft", c("All Diseased", "Incident that week"), lty=c(1,1), col=c("black", "red"))

#plot(x.vals, subpop.tot[11,]/subpop.grandtot, xlab="Year", ylab="Fraction of population", main="Death Rate", pch="")
#lines(x.vals, subpop.tot[11,]/subpop.grandtot)
#legend("topleft", c("Death Rate"), lty=c(1,1), col=c("black", "red"))


plot(x.vals, subpop.grandtot, ylim=c(0, max(subpop.grandtot)),  xlab="Year", ylab="People", main="Total Population", pch="")
lines(x.vals, subpop.grandtot)
dev.off()
}
options(scipen = 999)






