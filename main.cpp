#include <raylib.h>

#include <algorithm>
#include <cmath>
#include <vector>

int get_index(int x, int y, int n) { return x + y * n; }

class Fluid_cell
{
  int size;
  float dt;
  float diff;
  float visc;
  std::vector<float> s;
  std::vector<float> density;
  std::vector<float> Vx;
  std::vector<float> Vy;
  std::vector<float> Vx0;
  std::vector<float> Vy0;

  void set_bnd(int b, float *x)
  {
    for (int i = 1; i < size - 1; i++)
    {
      x[get_index(i, 0, size)] = b == 2 ? -x[get_index(i, 1, size)] : x[get_index(i, 1, size)];
      x[get_index(i, size - 1, size)] = b == 2 ? -x[get_index(i, size - 2, size)] : x[get_index(i, size - 2, size)];
    }
    for (int j = 1; j < size - 1; j++)
    {
      x[get_index(0, j, size)] = b == 1 ? -x[get_index(1, j, size)] : x[get_index(1, j, size)];
      x[get_index(size - 1, j, size)] = b == 1 ? -x[get_index(size - 2, j, size)] : x[get_index(size - 2, j, size)];
    }
    x[get_index(0, 0, size)] = 0.5f * (x[get_index(1, 0, size)] + x[get_index(0, 1, size)]);
    x[get_index(0, size - 1, size)] = 0.5f * (x[get_index(1, size - 1, size)] + x[get_index(0, size - 2, size)]);
    x[get_index(size - 1, 0, size)] = 0.5f * (x[get_index(size - 2, 0, size)] + x[get_index(size - 1, 1, size)]);
    x[get_index(size - 1, size - 1, size)] = 0.5f * (x[get_index(size - 2, size - 1, size)] + x[get_index(size - 1, size - 2, size)]);
  }

  void lin_solve(int b, float *x, float *x0, float a, float c)
  {
    float cRecip = 1.0f / c;
    for (int k = 0; k < 20; k++)
    {
      for (int j = 1; j < size - 1; j++)
      {
        for (int i = 1; i < size - 1; i++)
        {
          x[get_index(i, j, size)] = (x0[get_index(i, j, size)] + a * (x[get_index(i + 1, j, size)] + x[get_index(i - 1, j, size)] +
                                                                       x[get_index(i, j + 1, size)] + x[get_index(i, j - 1, size)])) *
                                     cRecip;
        }
      }
      set_bnd(b, x);
    }
  }

  void diffuse(int b, float *x, float *x0)
  {
    float a = dt * diff * (size - 2) * (size - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a);
  }

  void diffuse_vel(int b, float *x, float *x0)
  {
    float a = dt * visc * (size - 2) * (size - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a);
  }
  void project()
  {
    std::vector<float> div(size * size);
    std::vector<float> p(size * size);

    for (int j = 1; j < size - 1; j++)
    {
      for (int i = 1; i < size - 1; i++)
      {
        div[get_index(i, j, size)] = -0.5f *
                                     (Vx[get_index(i + 1, j, size)] - Vx[get_index(i - 1, j, size)] + Vy[get_index(i, j + 1, size)] -
                                      Vy[get_index(i, j - 1, size)]) /
                                     size;
        p[get_index(i, j, size)] = 0;
      }
    }

    set_bnd(0, div.data());
    set_bnd(0, p.data());
    lin_solve(0, p.data(), div.data(), 1, 6);

    for (int j = 1; j < size - 1; j++)
    {
      for (int i = 1; i < size - 1; i++)
      {
        Vx[get_index(i, j, size)] -= 0.5f * (p[get_index(i + 1, j, size)] - p[get_index(i - 1, j, size)]) * size;
        Vy[get_index(i, j, size)] -= 0.5f * (p[get_index(i, j + 1, size)] - p[get_index(i, j - 1, size)]) * size;
      }
    }
    set_bnd(1, Vx.data());
    set_bnd(2, Vy.data());
  }

  void compute_curl()
  {
    std::vector<float> curl(size * size);

    for (int j = 1; j < size - 1; j++)
    {
      for (int i = 1; i < size - 1; i++)
      {
        float dVxDy = (Vx[get_index(i, j + 1, size)] - Vx[get_index(i, j - 1, size)]) / (2.0f);  // Removed dt
        float dVyDx = (Vy[get_index(i + 1, j, size)] - Vy[get_index(i - 1, j, size)]) / (2.0f);  // Removed dt

        curl[get_index(i, j, size)] = dVyDx - dVxDy;
      }
    }

    // Calculate average velocity magnitude
    float avg_velocity_magnitude = 0.0f;
    for (int i = 0; i < size * size; ++i)
      avg_velocity_magnitude += std::sqrt(Vx[i] * Vx[i] + Vy[i] * Vy[i]);
    avg_velocity_magnitude /= (size * size);

    // Scale curl strength based on average velocity magnitude
    float max_curl_strength = 2.0f;  // Maximum curl strength
    float min_curl_strength = 0.2f;  // Minimum curl strength
    float curl_strength = min_curl_strength + (max_curl_strength - min_curl_strength) * (avg_velocity_magnitude);
    curl_strength = std::min(max_curl_strength, std::max(min_curl_strength, curl_strength));

    for (int j = 1; j < size - 1; j++)
    {
      for (int i = 1; i < size - 1; i++)
      {
        float dCurlX = (std::abs(curl[get_index(i + 1, j, size)]) - std::abs(curl[get_index(i - 1, j, size)])) / 2.0f;
        float dCurlY = (std::abs(curl[get_index(i, j + 1, size)]) - std::abs(curl[get_index(i, j - 1, size)])) / 2.0f;

        float length = std::sqrt(dCurlX * dCurlX + dCurlY * dCurlY) + 0.000001f;
        dCurlX /= length;
        dCurlY /= length;

        float curlValue = curl[get_index(i, j, size)];

        Vx[get_index(i, j, size)] += curl_strength * dCurlY * curlValue;
        Vy[get_index(i, j, size)] -= curl_strength * dCurlX * curlValue;
      }
    }
  }
  void advect(int b, float *d, float *d0, float *velocX, float *velocY)
  {
    float dtx = dt * (size - 2);
    float dty = dt * (size - 2);

    for (int j = 1; j < size - 1; j++)
    {
      for (int i = 1; i < size - 1; i++)
      {
        float x = i - dtx * velocX[get_index(i, j, size)];
        float y = j - dty * velocY[get_index(i, j, size)];

        x = std::max(0.5f, std::min(x, size - 1.5f));
        y = std::max(0.5f, std::min(y, size - 1.5f));

        int i0 = std::floor(x);
        int j0 = std::floor(y);
        int i1 = i0 + 1;
        int j1 = j0 + 1;

        float s1 = x - i0;
        float s0 = 1.0f - s1;
        float t1 = y - j0;
        float t0 = 1.0f - t1;

        d[get_index(i, j, size)] = s0 * (t0 * d0[get_index(i0, j0, size)] + t1 * d0[get_index(i0, j1, size)]) +
                                   s1 * (t0 * d0[get_index(i1, j0, size)] + t1 * d0[get_index(i1, j1, size)]);
      }
    }
    set_bnd(b, d);
  }

public:
  Fluid_cell(int size, float diffusion, float viscosity, float timestep) : size(size), diff(diffusion), visc(viscosity), dt(timestep)
  {
    s.resize(size * size);
    density.resize(size * size);
    Vx.resize(size * size);
    Vy.resize(size * size);
    Vx0.resize(size * size);
    Vy0.resize(size * size);
  }

  void step()
  {
    diffuse_vel(1, Vx0.data(), Vx.data());
    diffuse_vel(2, Vy0.data(), Vy.data());

    project();

    advect(1, Vx.data(), Vx0.data(), Vx.data(), Vy.data());
    advect(2, Vy.data(), Vy0.data(), Vx.data(), Vy.data());
    compute_curl();
    project();

    diffuse(0, s.data(), density.data());
    advect(0, density.data(), s.data(), Vx.data(), Vy.data());
  }

  void add_density(int x, int y, float amount)
  {
    int index = get_index(x, y, size);
    if (index >= 0 && index < size * size)
      density[index] += amount;
  }

  void add_velocity(int x, int y, float amountX, float amountY)
  {
    int index = get_index(x, y, size);
    if (index >= 0 && index < size * size)
    {
      Vx[index] += amountX;
      Vy[index] += amountY;
    }
  }

  void render(int screenWidth, int screenHeight)
  {
    int cellWidth = screenWidth / size;
    int cellHeight = screenHeight / size;

    for (int j = 0; j < size; j++)
    {
      for (int i = 0; i < size; i++)
      {
        float d = density[get_index(i, j, size)];

        float brightness_g = std::min(d, 1.0f) * 0.8f;
        float brightness_b = std::min(d, 1.0f) * 0.8f;
        Color cellColor = Color{0, static_cast<unsigned char>(brightness_g * 255), static_cast<unsigned char>(brightness_b * 255), 255};

        DrawRectangle(i * cellWidth, j * cellHeight, cellWidth, cellHeight, cellColor);
      }
    }
  }
};

int main()
{
  const int screenWidth = 800;
  const int screenHeight = 600;
  const int gridSize = 100;
  const float density = 16;
  const float vel_factor = 1;

  InitWindow(screenWidth, screenHeight, "Fluid Simulation - Enhanced");
  SetTargetFPS(60);
  Fluid_cell fluid(gridSize, 0.000003f, 0.0000007f, 0.2f);

  while (!WindowShouldClose())
  {
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT))
    {
      Vector2 mousePos = GetMousePosition();
      int x = static_cast<int>(mousePos.x / (screenWidth / gridSize));
      int y = static_cast<int>(mousePos.y / (screenHeight / gridSize));

      fluid.add_density(x, y, density);
      fluid.add_velocity(x, y, GetMouseDelta().x * vel_factor, GetMouseDelta().y * vel_factor);
    }

    fluid.step();

    BeginDrawing();
    ClearBackground(BLACK);

    fluid.render(screenWidth, screenHeight);

    EndDrawing();
  }

  CloseWindow();
  return 0;
}
