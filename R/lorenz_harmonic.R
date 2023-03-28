library(deSolve)
library(dplyr)

#
# Lorenz model coupled to harmonic oscillator.
#
lorenz_harmonic <- function(n = 1000,
                            skip = 100,
                            dt = 0.01,
                            coupling = 0.01) {
  # Internal function to compute derivatives for the ode solver.
  # x, y, z: The Lorenz system
  # u, v: position and velocity of harmonic oscillator, perturbed by Lorenz-x.
  # A force proportional to square of the x-component of the Lorenz model
  # is exerted on the oscillator with a coupling constant c.
  coupled <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dx <-  sigma * (y - x)
      dy <- x * (rho - z) - y
      dz <- x * y - beta * z
      du <- v
      dv <- -k * u + c * x**2
      return(list(c(dx, dy, dz, du, dv)))
    })
  }
  # Parameters for the model. The constant c is taken from the function
  # call to lorenz_harmonic, so it can be varied.
  params <- c(sigma = 10,
              rho = 28,
              beta = 8/3,
              k = 10,
              c = coupling)
  # Initial value of the variables.
  state_init <- c(x = 10, y = 10, z = 10, u = 10, v = 0)
  # A vector of time values created from the step-size, dt, the
  # number of points to skip (to avoid transient behavior) and
  # the number of points to output, n.
  times <- seq(0, dt * (skip + n), dt)
  # Use ode to solve the system
  ode_data <- ode(y = state_init, times = times,
                  func = coupled, parms = params)
  # Convert the data to a data frame
  model_data <- as.data.frame(ode_data)
  if (skip > 0) {
    model_data <- model_data[(skip + 1):(skip + n), ]
  }
  return(model_data)
}

# Function to extract the Lorenz part of the data as a matrix
# The input model_df is the output from lorenz_harmonic()
extract_lorenz <- function(model_df) {
  return(
    as.matrix(
      cbind(
        model_df |> pull(x),
        model_df |> pull(y),
        model_df |> pull(z)
      )
    )
  )
}

# Function to extract the oscillator part of the data as a matrix
# The input model_df is the output from lorenz_harmonic()
extract_oscillator <- function(model_df) {
  return(
    as.matrix(
      cbind(
        model_df |> pull(u),
        model_df |> pull(v)
      )
    )
  )
}
