#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592
#define L 200
#define R_Max 49
#define Adherent_P 0.95
#define Particles_No 1000

int main()
{
  int i, j, x, y, site[L][L], IsAdhered;
  float randAngle, randWalk, randAdherent;
  FILE *outfile;

  srand(time(0));
  outfile = fopen("out0.95.txt", "w");

  // --- site initialisation
  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++)
      site[i][j] = 0;
  site[L / 2][L / 2] = 1;

  // --- just letting the particles diffuse and adhere
  for (i = 0; i < Particles_No;)
  {
    // --- x and y initialisation
    randAngle = 2 * PI * ((float)rand()) / ((float)RAND_MAX);
    x = L / 2 + (int)(R_Max * cos(randAngle));
    y = L / 2 + (int)(R_Max * sin(randAngle));

    while (sqrt((x - L / 2) * (x - L / 2) + (y - L / 2) * (y - L / 2)) <= 2 * R_Max)
    {
      // --- checking the neighbouring sites
      if (site[x - 1][y] != 0 || site[x + 1][y] != 0 || site[x][y - 1] != 0 || site[x][y + 1] != 0)
      {
        randAdherent = ((float)rand()) / ((float)RAND_MAX);
        if (randAdherent <= Adherent_P)
        {
          site[x][y] = 1;
          i++;
          break;
        }
      }

      // --- random walking
      randWalk = ((float)rand()) / ((float)RAND_MAX);

      if (randWalk < 0.25) // Moving right
        if (site[x + 1][y] != 0)
          x--;
        else
          x++;
      else if (randWalk >= 0.25 && randWalk < 0.5) // Moving up
        if (site[x][y + 1] != 0)
          y--;
        else
          y++;
      else if (randWalk >= 0.5 && randWalk < 0.75) // Moving left
        if (site[x - 1][y] != 0)
          x++;
        else
          x--;
      else if (randWalk >= 0.75 && randWalk <= 1) // Moving down
        if (site[x][y - 1] != 0)
          y++;
        else
          y--;
    }
  }

  // --- printing the result to file
  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++)
      if (site[i][j] == 1)
        fprintf(outfile, "%d\t%d\n", i, j);

  return 0;
}
