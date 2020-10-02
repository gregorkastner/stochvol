## Simulate a highly persistent SV process of length 500
sim <- svsim(500, phi = 0.99, sigma = 0.1)

print(sim)
summary(sim)
plot(sim)

## Simulate an SV process with leverage
sim <- svsim(200, phi = 0.94, sigma = 0.15, rho = -0.6)

print(sim)
summary(sim)
plot(sim)

## Simulate an SV process with conditionally heavy-tails
sim <- svsim(250, phi = 0.91, sigma = 0.05, nu = 5)

print(sim)
summary(sim)
plot(sim)

