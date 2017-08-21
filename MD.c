#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 40 // Number of particles
#define L 10 // Length of the box
#define equilT 10 // Equilibrium temprature

double
    x[N],
    y[N],               // Coordinates of particles
    vx[N], vy[N],       // Velocity of particles
    x_ax[N], x_ay[N],   // Ex-Acceleration of particles
    ax[N], ay[N],       // Acceleration of particles
    fx[N], fy[N],       // Total exerted force on each particle
    rx[N][N], ry[N][N]; // The distance of particle j from particle i

double LJ_Force(double r);
double LJ_Potential(double r);
double calTemperature();
void calForces();
void calDistances();
void initialConditions();
void verlet(double dt, double T);
void calAccelerations();
void cold(double currentT, double futureT);

int main()
{
  verlet(0.001, 30 * equilT);

  return 0;
}

void initialConditions()
{
  double v0xCM, v0yCM; // Initial velocity of Center of Mass
  int i, j;

  // Setting initial coordinates
  for (i = 0; i < 5; i++)
  {
    for (j = 0; j < 8; j++)
    {
      x[8 * i + j] = 1 + j;
      y[8 * i + j] = 1 + 2 * i;
    }
  }
  /* Initial configuration of particles in the left side of the box
      - - - - - - - - - -
     |                   |
     |  * * * * * * * *  |
     |                   |
     |  * * * * * * * *  |
     |                   |
     |  * * * * * * * *  |
     |                   |
     |  * * * * * * * *  |
     |                   |
     |  * * * * * * * *  |
     |                   |
      - - - - - - - - - -   */

  // Setting particles initial velocities to random numbers and CM initial velocity to zero
  v0xCM = 0;
  v0yCM = 0;

  for (i = 0; i < N; i++)
  {
    vx[i] = ((double)(rand() % 1000)) / 999;
    vy[i] = ((double)(rand() % 1000)) / 999;

    v0xCM += vx[i];
    v0yCM += vy[i];
  }

  v0xCM = v0xCM / N;
  v0yCM = v0yCM / N;

  for (i = 0; i < N; i++)
  {
    vx[i] -= v0xCM;
    vy[i] -= v0yCM;
  }
}

void calDistances()
{
  int i, j;

  // Calculating the distances
  for (i = 0; i < N; i++)
  {
    for (j = i + 1; j < N; j++)
    {
      rx[i][j] = x[j] - x[i];
      if (rx[i][j] > (L / 2))
        rx[i][j] -= L;
      if (rx[i][j] < ((-1) * L / 2))
        rx[i][j] += L;

      ry[i][j] = y[j] - y[i];
      if (ry[i][j] > (L / 2))
        ry[i][j] -= L;
      if (ry[i][j] < ((-1) * L / 2))
        ry[i][j] += L;
    }
  }
}

void calForces()
{
  int i, j;
  double r_temp, fx_temp, fy_temp;

  // Calculating forces exerted on each particle
  for (i = 0; i < N; i++)
  {
    fx[i] = 0;
    fy[i] = 0;
  }

  for (i = 0; i < N; i++)
  {
    for (j = i + 1; j < N; j++)
    {
      r_temp = sqrt(pow(rx[i][j], 2) + pow(ry[i][j], 2));
      fx_temp = LJ_Force(r_temp) * rx[i][j] / r_temp;
      fx[i] -= fx_temp;
      fx[j] += fx_temp;
      fy_temp = LJ_Force(r_temp) * ry[i][j] / r_temp;
      fy[i] -= fy_temp;
      fy[j] += fy_temp;
    }
  }
}

void calAccelerations()
{
  int i;

  for (i = 0; i < N; i++)
  {
    x_ax[i] = ax[i];
    x_ay[i] = ay[i];
    ax[i] = fx[i];
    ay[i] = fy[i];
  }
}

void verlet(double dt, double T)
{
  int i, j;
  double t = 0, temperature, tempT = 0.8;
  FILE *outfile;

  initialConditions();

  calDistances();
  calForces();
  calAccelerations();

  outfile = fopen("Temperatures.txt", "w");
  fprintf(outfile, "%f\t%f\n", t, calTemperature());

  for (i = 0; i < ((int)(T / dt)); i++)
  {
    t += dt;

    for (j = 0; j < N; j++)
    {
      x[j] += vx[j] * dt + 0.5 * ax[j] * pow(dt, 2);
      if (x[j] > L)
        x[j] -= L;
      if (x[j] < 0)
        x[j] += L;
      y[j] += vy[j] * dt + 0.5 * ay[j] * pow(dt, 2);
      if (y[j] > L)
        y[j] -= L;
      if (y[j] < 0)
        y[j] += L;
    }

    calDistances();
    calForces();
    calAccelerations();

    for (j = 0; j < N; j++)
    {
      vx[j] += 0.5 * (x_ax[j] + ax[j]) * dt;
      vy[j] += 0.5 * (x_ay[j] + ay[j]) * dt;
    }

    temperature = calTemperature();
    fprintf(outfile, "%f\t%f\n", t, calTemperature());

    if ((i % ((int)(equilT / dt))) == 0)
    {
      if (i > 0 && i <= (3 * ((int)(equilT / dt))))
        cold(temperature, 0.8);
      else
      {
        tempT = (tempT - 0.005) / 2 + 0.005;
        cold(temperature, tempT);
      }
    }
  }

  fclose(outfile);

  FILE *places, *velocities;
  places = fopen("finalPlaces.txt", "w");
  velocities = fopen("finalVelocities.txt", "w");

  for (i = 0; i < N; i++)
  {
    fprintf(places, "%f\t%f\n", x[i], y[i]);
    fprintf(velocities, "%f\t%f\n", vx[i], vy[i]);
  }

  fclose(places);
  fclose(velocities);
}

void cold(double currentT, double futureT)
{
  int i = 0;

  for (i = 0; i < N; i++)
  {
    vx[i] *= sqrt(futureT / currentT);
    vy[i] *= sqrt(futureT / currentT);
  }
}

double calTemperature()
{
  double vxSquaredMean = 0, vySquaredMean = 0, temperature;
  int i;

  for (i = 0; i < N; i++)
  {
    vxSquaredMean += pow(vx[i], 2);
    vySquaredMean += pow(vy[i], 2);
  }

  vxSquaredMean /= N;
  vySquaredMean /= N;

  temperature = 0.5 * (vxSquaredMean + vySquaredMean) / (N - 1);

  return temperature;
}

double LJ_Potential(double r)
{
  if (r <= 2)
    return (4 / pow(r, 12) - 4 / pow(r, 6) - 0.181640625 * r + 0.4248046875);
  else
    return 0;
}

double LJ_Force(double r)
{
  if (r <= 2)
    return (48 / pow(r, 13) - 24 / pow(r, 7) + 0.181640625);
  else
    return 0;
}
