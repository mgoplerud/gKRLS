library(gKRLS)

# Should RUN
N <- 100
x1 <- rnorm(N)
x2 <- rexp(N)
kern_data <- cbind(x1 , x2)
fake_id <- sample(id, size = N, replace = T)
fake_outcome <- rnorm(N, sd = 10) + rnorm(26)[match(fake_id, letters)]
reg_data <- data.frame(outcome_y = fake_outcome, id = fake_id)
gKRLS_random <- gKRLS(formula = outcome_y ~ 1 + (1| id), 
                      kernel_X = kern_data, sketch_method = 'none',
                      data = reg_data, 
                      family = gaussian())


reg_data <- data.frame(response = fake_outcome, id = fake_id)
gKRLS_random2 <- gKRLS(formula = response ~ 1 + (1| id), 
                      kernel_X = kern_data, 
                      data = reg_data, sketch_method = 'none',
                      family = gaussian())


gKRLS_random3 <- gKRLS(formula = fake_outcome ~ 1 + (1| fake_id), 
                      kernel_X = kern_data, sketch_method = 'none',
                      data = NULL, 
                      family = gaussian())

all.equal(gKRLS_random3, gKRLS_random)

