// EPSILON is a tolerance value for floating point roundoff error.
// It is used in many calculations where we want to err
// on a certain side of a threshold, such as determining
// whether or not a point is inside a solid or not,
// or whether a point is at least a minimum distance
// away from another point.
const float kEPSILON = 1.0e-6f;