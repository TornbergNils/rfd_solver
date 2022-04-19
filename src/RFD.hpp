#ifndef RFD_H
#define RFD_H

#include <cmath>
#include <vector>

// Dimension indicates what part of cross prod, 0,1,2 => x,y,z
double Get_cross_product(const int ix, const EM_field_matrix &EM_field,
                         const int dimension)
{
  double cross_product;

  // x-component of cross prod
  if (dimension == 0)
  {

    cross_product = (EM_field.E_y[ix] * EM_field.B_z[ix] -
                     EM_field.E_z[ix] * EM_field.B_y[ix]);
    // y-component
  }
  else if (dimension == 1)
  {
    cross_product = -(EM_field.E_x[ix] * EM_field.B_z[ix] -
                      EM_field.E_z[ix] * EM_field.B_x[ix]);
    // z-component
  }
  else if (dimension == 2)
  {
    cross_product = (EM_field.E_x[ix] * EM_field.B_y[ix] -
                     EM_field.E_y[ix] * EM_field.B_x[ix]);
  }
  else
  {
    printf("Invalid dimension for cross product!");
    cross_product = 0.0;
  }

  return cross_product;
}

double get_w(const double E_cross_B_squared, const double E_squared,
             const double B_squared)
{

  double E_sq_plus_B_sq = E_squared + B_squared;
  double w = 4.0 * E_cross_B_squared / (E_sq_plus_B_sq * E_sq_plus_B_sq);
  return w;
}

double get_u(const double w, const double E_squared, const double B_squared,
             const double eps)
{

  double factor1 = 2.0 * B_squared / (E_squared + B_squared);
  double factor2;
  if (w > 1.0)
  {
    factor2 = 1.0;
  }
  else if (w > eps)
  {
    factor2 = (1.0 - std::sqrt(1 - w)) / w;
  }
  else
  {
    factor2 = 1.0 / 2.0 + w / 8;
  }
  double u = factor1 * factor2;
  // Return u only if less than 1, otherwise 1
  return u * ( u < 1) + ( u > 1);
}

double get_RFD_component(const double u, const double w,
                         const double E_cross_B_component,
                         const double B_component, const double E_component,
                         const double E_squared, const double B_squared,
                         const double E_dot_B, const double E_cross_B_squared,
                         const double sign)
{

  const double term1 = std::sqrt(u - u * u) * E_cross_B_component;

  const double term2 = (1 - u) * std::sqrt(B_squared) * E_component +
                       u * E_dot_B * B_component / std::sqrt(B_squared);

  const double factor1 = term1 + sign * term2;

  // to normalize analytically
  // const double factor2 =
  //    1.0 / std::sqrt(E_squared * B_squared - u * E_cross_B_squared);

  const double RFD_component = factor1; //* factor2 ;
  return RFD_component;
}

class RFD_matrix
{
public:
  int nx;
  int ny;
  std::vector<double> RFD_x;
  std::vector<double> RFD_y;
  std::vector<double> RFD_z;
  int sign;

  void Save(std::string filename, bool append)
  {
    std::ofstream savestream_x;
    std::ofstream savestream_y;
    std::ofstream savestream_z;
    if (append == true)
    {
      savestream_x.open(filename + "_x",
                        std::ios::out | std::ios::app | std::ios::binary);
      savestream_y.open(filename + "_y",
                        std::ios::out | std::ios::app | std::ios::binary);
      savestream_z.open(filename + "_z",
                        std::ios::out | std::ios::app | std::ios::binary);
    }
    else
    {
      savestream_x.open(filename + "_x",
                        std::ios::out | std::ios::trunc | std::ios::binary);
      savestream_y.open(filename + "_y",
                        std::ios::out | std::ios::trunc | std::ios::binary);
      savestream_z.open(filename + "_z",
                        std::ios::out | std::ios::trunc | std::ios::binary);
    }

    savestream_x.write((char *)&RFD_x[0], RFD_x.size() * sizeof(double));
    savestream_y.write((char *)&RFD_y[0], RFD_y.size() * sizeof(double));
    savestream_z.write((char *)&RFD_z[0], RFD_z.size() * sizeof(double));
    savestream_x.close();
    savestream_y.close();
    savestream_z.close();
  }

  RFD_matrix(EM_field_matrix EM_field, int sign_init)
  {
    nx = EM_field.nx;
    ny = EM_field.ny;
    sign = sign_init;
    std::vector<double> temp_RFD_x(nx * ny);
    std::vector<double> temp_RFD_y(nx * ny);
    std::vector<double> temp_RFD_z(nx * ny);

    const double eps = 10e-15;

    for (int ix = 0; ix < nx * ny; ix++)
    {
      const double E_cross_B_x = Get_cross_product(ix, EM_field, 0);
      const double E_cross_B_y = Get_cross_product(ix, EM_field, 1);
      const double E_cross_B_z = Get_cross_product(ix, EM_field, 2);

      const double E_squared = EM_field.E_x[ix] * EM_field.E_x[ix] +
                               EM_field.E_y[ix] * EM_field.E_y[ix] +
                               EM_field.E_z[ix] * EM_field.E_z[ix];

      const double B_squared = EM_field.B_x[ix] * EM_field.B_x[ix] +
                               EM_field.B_y[ix] * EM_field.B_y[ix] +
                               EM_field.B_z[ix] * EM_field.B_z[ix];

      const double E_cross_B_squared = E_cross_B_x * E_cross_B_x +
                                       E_cross_B_y * E_cross_B_y +
                                       E_cross_B_z * E_cross_B_z;

      const double E_dot_B = EM_field.E_x[ix] * EM_field.B_x[ix] +
                             EM_field.E_y[ix] * EM_field.B_y[ix] +
                             EM_field.E_z[ix] * EM_field.B_z[ix];

      // this is obviusly unreasonable but avoid division by zero
      if ((E_squared + B_squared) < eps)
      {
        temp_RFD_x[ix] = 0.0;
        temp_RFD_y[ix] = 0.0;
        temp_RFD_z[ix] = 0.0;
      }
      // This case, E=0 needs more careful investigation
      else if (E_squared < eps)
      {
        // use copysign as signum

        temp_RFD_x[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.B_x[ix];
        temp_RFD_y[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.B_y[ix];
        temp_RFD_z[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.B_z[ix];

        double denominator = std::sqrt(temp_RFD_x[ix] * temp_RFD_x[ix] +
                                       temp_RFD_y[ix] * temp_RFD_y[ix] +
                                       temp_RFD_z[ix] * temp_RFD_z[ix]);
        temp_RFD_x[ix] /= denominator;
        temp_RFD_y[ix] /= denominator;
        temp_RFD_z[ix] /= denominator;
        // This case, B=0 needs more careful investigation
      }
      else if (B_squared < eps)
      {

        temp_RFD_x[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.E_x[ix];
        temp_RFD_y[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.E_y[ix];
        temp_RFD_z[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.E_z[ix];

        double denominator = std::sqrt(temp_RFD_x[ix] * temp_RFD_x[ix] +
                                       temp_RFD_y[ix] * temp_RFD_y[ix] +
                                       temp_RFD_z[ix] * temp_RFD_z[ix]);
        temp_RFD_x[ix] /= denominator;
        temp_RFD_y[ix] /= denominator;
        temp_RFD_z[ix] /= denominator;
        // Case 2b, E_dot_B = 0
      }
      else if ((E_dot_B * E_dot_B) < eps && E_squared <= B_squared + eps)
      {
        // printf( "%lf\n", E_dot_B );
        // Case 2bi, E_dot_B = 0 and |E| = |B|
        if (std::abs(std::sqrt(E_squared) - std::sqrt(B_squared)) < eps)
        {

          temp_RFD_x[ix] = E_cross_B_x;
          temp_RFD_y[ix] = E_cross_B_y;
          temp_RFD_z[ix] = E_cross_B_z;

          double denominator = std::sqrt(temp_RFD_x[ix] * temp_RFD_x[ix] +
                                         temp_RFD_y[ix] * temp_RFD_y[ix] +
                                         temp_RFD_z[ix] * temp_RFD_z[ix]);
          temp_RFD_x[ix] /= denominator;
          temp_RFD_y[ix] /= denominator;
          temp_RFD_z[ix] /= denominator;
          // Case 2bii, E_dot_B = 0 and E^2 < B^2
        }
        else
        {

          const double w = get_w(E_cross_B_squared, E_squared, B_squared);
          const double u = get_u(w, E_squared, B_squared, eps);

          double term1x = E_cross_B_x;
          double term1y = E_cross_B_y;
          double term1z = E_cross_B_z;

          double term2x;
          double term2y;
          double term2z;
          if (u > 1.0 - eps)
          {
            term2x = 0.0;
            term2y = 0.0;
            term2z = 0.0;
          }
          else
          {
            term2x = sign * std::sqrt(1.0 - u) * std::sqrt(B_squared) *
                     EM_field.E_x[ix];
            term2y = sign * std::sqrt(1.0 - u) * std::sqrt(B_squared) *
                     EM_field.E_y[ix];
            term2z = sign * std::sqrt(1.0 - u) * std::sqrt(B_squared) *
                     EM_field.E_z[ix];
          }
          double term3x = sign * EM_field.B_x[ix];
          double term3y = sign * EM_field.B_y[ix];
          double term3z = sign * EM_field.B_z[ix];

          double numeratorx = term1x + term2x + term3x;
          double numeratory = term1y + term2y + term3y;
          double numeratorz = term1z + term2z + term3z;

          temp_RFD_x[ix] = numeratorx;
          temp_RFD_y[ix] = numeratory;
          temp_RFD_z[ix] = numeratorz;

          double denominator = std::sqrt(temp_RFD_x[ix] * temp_RFD_x[ix] +
                                         temp_RFD_y[ix] * temp_RFD_y[ix] +
                                         temp_RFD_z[ix] * temp_RFD_z[ix]);
          temp_RFD_x[ix] /= denominator;
          temp_RFD_y[ix] /= denominator;
          temp_RFD_z[ix] /= denominator;
        }
      }
      else
      {
        const double w = get_w(E_cross_B_squared, E_squared, B_squared);

        const double u = get_u(w, E_squared, B_squared, eps);

        temp_RFD_x[ix] = get_RFD_component(
            u, w, E_cross_B_x, EM_field.B_x[ix], EM_field.E_x[ix], E_squared,
            B_squared, E_dot_B, E_cross_B_squared, sign);

        temp_RFD_y[ix] = get_RFD_component(
            u, w, E_cross_B_y, EM_field.B_y[ix], EM_field.E_y[ix], E_squared,
            B_squared, E_dot_B, E_cross_B_squared, sign);

        temp_RFD_z[ix] = get_RFD_component(
            u, w, E_cross_B_z, EM_field.B_z[ix], EM_field.E_z[ix], E_squared,
            B_squared, E_dot_B, E_cross_B_squared, sign);

        double denominator = std::sqrt(temp_RFD_x[ix] * temp_RFD_x[ix] +
                                       temp_RFD_y[ix] * temp_RFD_y[ix] +
                                       temp_RFD_z[ix] * temp_RFD_z[ix]);
        temp_RFD_x[ix] /= denominator;
        temp_RFD_y[ix] /= denominator;
        temp_RFD_z[ix] /= denominator;
      }
         if( std::isnan(temp_RFD_x[ix]) 
         || std::isnan(temp_RFD_y[ix]) 
         || std::isnan(temp_RFD_z[ix]) ) {printf("RFD nan at %d", ix); }
    }

    RFD_x = temp_RFD_x;
    RFD_y = temp_RFD_y;
    RFD_z = temp_RFD_z;
  }

  void Update(EM_field_matrix EM_field, int sign_init)
  {
    nx = EM_field.nx;
    ny = EM_field.ny;
    sign = sign_init;

    const double eps = 10e-15;

    for (int ix = 0; ix < nx * ny; ix++)
    {
      const double E_cross_B_x = Get_cross_product(ix, EM_field, 0);
      const double E_cross_B_y = Get_cross_product(ix, EM_field, 1);
      const double E_cross_B_z = Get_cross_product(ix, EM_field, 2);

      const double E_squared = EM_field.E_x[ix] * EM_field.E_x[ix] +
                               EM_field.E_y[ix] * EM_field.E_y[ix] +
                               EM_field.E_z[ix] * EM_field.E_z[ix];

      const double B_squared = EM_field.B_x[ix] * EM_field.B_x[ix] +
                               EM_field.B_y[ix] * EM_field.B_y[ix] +
                               EM_field.B_z[ix] * EM_field.B_z[ix];

      const double E_cross_B_squared = E_cross_B_x * E_cross_B_x +
                                       E_cross_B_y * E_cross_B_y +
                                       E_cross_B_z * E_cross_B_z;

      const double E_dot_B = EM_field.E_x[ix] * EM_field.B_x[ix] +
                             EM_field.E_y[ix] * EM_field.B_y[ix] +
                             EM_field.E_z[ix] * EM_field.B_z[ix];

      // this is obviusly unreasonable but avoid division by zero
      if ((E_squared + B_squared) < eps)
      {
        RFD_x[ix] = 0.0;
        RFD_y[ix] = 0.0;
        RFD_z[ix] = 0.0;
      }
      // This case, E=0 needs more careful investigation
      if (E_squared < eps)
      {
        // use copysign as signum

        RFD_x[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.B_x[ix];
        RFD_y[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.B_y[ix];
        RFD_z[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.B_z[ix];
      double denominator =
          std::sqrt(RFD_x[ix] * RFD_x[ix] + RFD_y[ix] * RFD_y[ix] +
                    RFD_z[ix] * RFD_z[ix]);
      RFD_x[ix] /= denominator;
      RFD_y[ix] /= denominator;
      RFD_z[ix] /= denominator;

        // This case, B=0 needs more careful investigation
      }
      else if (B_squared < eps)
      {

        RFD_x[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.E_x[ix];
        RFD_y[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.E_y[ix];
        RFD_z[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.E_z[ix];
      double denominator =
          std::sqrt(RFD_x[ix] * RFD_x[ix] + RFD_y[ix] * RFD_y[ix] +
                    RFD_z[ix] * RFD_z[ix]);
      RFD_x[ix] /= denominator;
      RFD_y[ix] /= denominator;
      RFD_z[ix] /= denominator;

        // Case 2b, E_dot_B = 0
      }
      else if ((E_dot_B * E_dot_B) < eps && E_squared <= B_squared + eps)
      {
        // printf( "%lf\n", E_dot_B );
        // Case 2bi, E_dot_B = 0 and |E| = |B|
        if (std::abs(std::sqrt(E_squared) - std::sqrt(B_squared)) < eps)
        {

          RFD_x[ix] = E_cross_B_x;
          RFD_y[ix] = E_cross_B_y;
          RFD_z[ix] = E_cross_B_z;
      double denominator =
          std::sqrt(RFD_x[ix] * RFD_x[ix] + RFD_y[ix] * RFD_y[ix] +
                    RFD_z[ix] * RFD_z[ix]);
      RFD_x[ix] /= denominator;
      RFD_y[ix] /= denominator;
      RFD_z[ix] /= denominator;
          // Case 2bii, E_dot_B = 0 and E^2 < B^2
        }
        else
        {

          const double w = get_w(E_cross_B_squared, E_squared, B_squared);
          const double u = get_u(w, E_squared, B_squared, eps);

          double term1x = E_cross_B_x;
          double term1y = E_cross_B_y;
          double term1z = E_cross_B_z;

          double term2x;
          double term2y;
          double term2z;
          if (u > 1.0 - eps)
          {
            term2x = 0.0;
            term2y = 0.0;
            term2z = 0.0;
          }
          else
          {
            term2x = sign * std::sqrt(1.0 - u) * std::sqrt(B_squared) *
                     EM_field.E_x[ix];
            term2y = sign * std::sqrt(1.0 - u) * std::sqrt(B_squared) *
                     EM_field.E_y[ix];
            term2z = sign * std::sqrt(1.0 - u) * std::sqrt(B_squared) *
                     EM_field.E_z[ix];
          }
          double term3x = sign * EM_field.B_x[ix];
          double term3y = sign * EM_field.B_y[ix];
          double term3z = sign * EM_field.B_z[ix];

          double numeratorx = term1x + term2x + term3x;
          double numeratory = term1y + term2y + term3y;
          double numeratorz = term1z + term2z + term3z;

          RFD_x[ix] = numeratorx;
          RFD_y[ix] = numeratory;
          RFD_z[ix] = numeratorz;
      double denominator =
          std::sqrt(RFD_x[ix] * RFD_x[ix] + RFD_y[ix] * RFD_y[ix] +
                    RFD_z[ix] * RFD_z[ix]);
      RFD_x[ix] /= denominator;
      RFD_y[ix] /= denominator;
      RFD_z[ix] /= denominator;
        }

      }
      else
      {
        const double w = get_w(E_cross_B_squared, E_squared, B_squared);

        const double u = get_u(w, E_squared, B_squared, eps);

        RFD_x[ix] = get_RFD_component(u, w, E_cross_B_x, EM_field.B_x[ix],
                                      EM_field.E_x[ix], E_squared, B_squared,
                                      E_dot_B, E_cross_B_squared, sign);

        RFD_y[ix] = get_RFD_component(u, w, E_cross_B_y, EM_field.B_y[ix],
                                      EM_field.E_y[ix], E_squared, B_squared,
                                      E_dot_B, E_cross_B_squared, sign);

        RFD_z[ix] = get_RFD_component(u, w, E_cross_B_z, EM_field.B_z[ix],
                                      EM_field.E_z[ix], E_squared, B_squared,
                                      E_dot_B, E_cross_B_squared, sign);
      double denominator =
          std::sqrt(RFD_x[ix] * RFD_x[ix] + RFD_y[ix] * RFD_y[ix] +
                    RFD_z[ix] * RFD_z[ix]);
      RFD_x[ix] /= denominator;
      RFD_y[ix] /= denominator;
      RFD_z[ix] /= denominator;
      }
    }
  }
};

#endif // RFD_H