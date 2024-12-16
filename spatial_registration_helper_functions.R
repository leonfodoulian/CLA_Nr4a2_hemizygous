# Calculate angle of a line using atan2 and line endpoint coordinates
CalculateLineAngle <- function(
    x1,
    y1,
    x2,
    y2
){
  return(atan2(y = y2 - y1,
               x = x2 - x1))
}

# Adjust angle by pi (flips angles to account for lines with opposite directions)
AdjustAngleByPI <- function(
    angle
){
  return(angle + pi)
}

# Normalize angle to [0, 2pi) range
NormalizeAngle02PI <- function(
    angle
){
  return(angle %% (2 * pi))
}

# Compute the angle in radian between two vectors
ComputeRotationAngle <- function(
    query.coords,
    reference.coords,
    x1.name,
    x2.name,
    y1.name,
    y2.name
){
  # Get angle of query data
  query.angle <- NormalizeAngle02PI(
    angle = CalculateLineAngle(
      x1 = query.coords[[x1.name]],
      y1 = query.coords[[y1.name]],
      x2 = query.coords[[x2.name]],
      y2 = query.coords[[y2.name]]
    )
  )
  
  # Get angle of reference data
  reference.angle <- NormalizeAngle02PI(
    angle = CalculateLineAngle(
      x1 = reference.coords[[x1.name]],
      y1 = reference.coords[[y1.name]],
      x2 = reference.coords[[x2.name]],
      y2 = reference.coords[[y2.name]]
    )
  )
  
  # Adjust query angle if lines are in opposite directions
  if (abs(x = reference.angle - query.angle) > (pi/2) && abs(x = reference.angle - query.angle) < (3*pi/2)) {
    query.angle <- NormalizeAngle02PI(
      angle = AdjustAngleByPI(
        angle = query.angle
      )
    )
  }
  
  # Get angle between the two vectors
  rotation.angle <- NormalizeAngle02PI(angle = reference.angle - query.angle)
  
  # Return rotation angle
  return(rotation.angle)
}

# Perform rotation of axes in two dimensions
RotateCoordinates <- function(
    xy,
    angle
){
  # Transform data to matrix
  xy <- as.matrix(x = xy)
  
  # Find cos and sin of the angle
  cos.angle <- cos(x = angle)
  sin.angle <- sin(x = angle)
  
  # Rotate coordinates
  xy.rot <- xy %*% t(x = matrix(data = c(cos.angle, sin.angle, -sin.angle, cos.angle),
                                nrow = 2,
                                ncol = 2))
  
  # Return rotated coordinates
  return(xy.rot)
}
