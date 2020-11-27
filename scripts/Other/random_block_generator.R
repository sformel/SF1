#random generator for blocks in Chapter 3

setwd("~/Google Drive/Van Bael Google Drive/Boggs443 Data & User folders/Users/Grad Students/Steve/Chapter3_seed_death/")

A <- 1:12
B <- 13:24
C <- 25:36
D <- 37:48

B1 <- c(A[1:3],B[1:3],C[1:3],D[1:3])
B2 <- c(A[4:6],B[4:6],C[4:6],D[4:6])
B3 <- c(A[7:9],B[7:9],C[7:9],D[7:9])
B4 <- c(A[10:12],B[10:12],C[10:12],D[10:12])

R1 <- sample(B1)
R2 <- sample(B2)
R3 <- sample(B3)
R4 <- sample(B4)

blocks <- data.frame("Block 1" = R1,"Block 2" =R2, "Block 3" =R3, "Block 4" =R4)

write.table(blocks, file = "chap3_randomized_blocks.txt")


#for greenhouse sampling
A <- 1:12
B <- 13:24
C <- 25:36
D <- 37:48

B1 <- c(A[1:3],B[1:3],C[1:3],D[1:3])
B2 <- c(A[4:6],B[4:6],C[4:6],D[4:6])
B3 <- c(A[7:9],B[7:9],C[7:9],D[7:9])
B4 <- c(A[10:12],B[10:12],C[10:12],D[10:12])

R1 <- sample(B1)
R2 <- sample(B2)
R3 <- sample(B3)
R4 <- sample(B4)

blocks <- data.frame("Block 1" = R1,"Block 2" =R2, "Block 3" =R3, "Block 4" =R4)

write.table(blocks, file = "chap3_randomized_blocks.txt")