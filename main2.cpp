#include <raylib.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

int get_index(int x, int y, int n) { return x + y * n; }

class FluidCube
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
  FluidCube(int size, float diffusion, float viscosity, float timestep) : size(size), diff(diffusion), visc(viscosity), dt(timestep)
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
    diffuse(1, Vx0.data(), Vx.data());
    diffuse(2, Vy0.data(), Vy.data());

    project();

    advect(1, Vx.data(), Vx0.data(), Vx.data(), Vy.data());
    advect(2, Vy.data(), Vy0.data(), Vx.data(), Vy.data());

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

  void render(int screenWidth, int screenHeight, bool showVelocity)
  {
    int cellWidth = screenWidth / size;
    int cellHeight = screenHeight / size;

    for (int j = 0; j < size; j++)
    {
      for (int i = 0; i < size; i++)
      {
        Color cellColor;
        if (showVelocity)
        {
          float vx = Vx[get_index(i, j, size)];
          float vy = Vy[get_index(i, j, size)];

          float res = std::sqrt(vx * vx + vy * vy);

          float brightness_g = std::min(res, 1.0f) * 0.8f;
          float brightness_b = std::min(res, 1.0f) * 0.5f;
          cellColor = Color{0, static_cast<unsigned char>(brightness_g * 255), static_cast<unsigned char>(brightness_b * 255), 255};
        }
        else
        {
          float d = density[get_index(i, j, size)];
          float brightness_g = std::min(d, 1.0f) * 0.8f;
          float brightness_b = std::min(d, 1.0f) * 0.5f;
          cellColor = Color{0, static_cast<unsigned char>(brightness_g * 255), static_cast<unsigned char>(brightness_b * 255), 255};
        }
        DrawRectangle(i * cellWidth, j * cellHeight, cellWidth, cellHeight, cellColor);
      }
    }
  }
  void set_diffusion(float diffusion) { diff = diffusion; }
  void set_viscosity(float viscosity) { visc = viscosity; }
  void set_timestep(float timestep) { dt = timestep; }
};

// Helper function to convert float to string with limited precision
std::string float_to_string(float value)
{
  std::ostringstream out;
  out << std::fixed << std::setprecision(6) << value;
  return out.str();
}

int main()
{
  const int screenWidth = 1200;  // Increased width to accommodate controls
  const int screenHeight = 600;
  const int gridSize = 100;
  const float baseDensity = 0.5f;
  const float baseVelFactor = 0.5f;

  InitWindow(screenWidth, screenHeight, "Fluid Simulation - Enhanced");
  SetTargetFPS(60);

  // Initial fluid parameters
  float diffusion = 0.000005f;
  float viscosity = 0.00001f;
  float timestep = 0.2f;
  float densityAmount = baseDensity;
  float velocityFactor = baseVelFactor;
  int brushRadius = 2;
  bool showVelocity = false;

  FluidCube fluid(gridSize, diffusion, viscosity, timestep);

  // GUI Parameters
  int selectedParam = 0;
  const char *paramNames[] = {"Diffusion", "Viscosity", "Timestep", "Density", "Velocity Factor", "Brush Radius"};
  const int paramCount = sizeof(paramNames) / sizeof(paramNames[0]);

  while (!WindowShouldClose())
  {
    // Parameter adjustment controls
    if (IsKeyPressed(KEY_UP))
      selectedParam = (selectedParam - 1 + paramCount) % paramCount;
    if (IsKeyPressed(KEY_DOWN))
      selectedParam = (selectedParam + 1) % paramCount;

    // Parameter value adjustment
    float adjustmentFactor = IsKeyDown(MOUSE_BUTTON_LEFT) ? 1.1f : 0.9f;
    switch (selectedParam)
    {
      case 0:  // Diffusion
        diffusion *= IsKeyPressed(KEY_RIGHT) ? adjustmentFactor : (IsKeyPressed(KEY_LEFT) ? 1.0f / adjustmentFactor : 1.0f);
        fluid.set_diffusion(diffusion);
        break;
      case 1:  // Viscosity
        viscosity *= IsKeyPressed(KEY_RIGHT) ? adjustmentFactor : (IsKeyPressed(KEY_LEFT) ? 1.0f / adjustmentFactor : 1.0f);
        fluid.set_viscosity(viscosity);
        break;
      case 2:  // Timestep
        timestep *= IsKeyPressed(KEY_RIGHT) ? adjustmentFactor : (IsKeyPressed(KEY_LEFT) ? 1.0f / adjustmentFactor : 1.0f);
        fluid.set_timestep(timestep);
        break;
      case 3:  // Density Amount
        densityAmount *= IsKeyPressed(KEY_RIGHT) ? adjustmentFactor : (IsKeyPressed(KEY_LEFT) ? 1.0f / adjustmentFactor : 1.0f);
        break;
      case 4:  // Velocity Factor
        velocityFactor *= IsKeyPressed(KEY_RIGHT) ? adjustmentFactor : (IsKeyPressed(KEY_LEFT) ? 1.0f / adjustmentFactor : 1.0f);
        break;
      case 5:  // Brush Radius
        brushRadius += IsKeyPressed(KEY_RIGHT) ? 1 : (IsKeyPressed(KEY_LEFT) ? -1 : 0);
        brushRadius = std::max(1, std::min(brushRadius, 10));
        break;
    }

    // Toggle velocity/density view
    if (IsKeyPressed(KEY_V))
      showVelocity = !showVelocity;

    // Fluid interaction
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT))
    {
      Vector2 mousePos = GetMousePosition();

      // Only interact if mouse is within the simulation area (excluding control panel)
      if (mousePos.x < screenWidth - 200)
      {
        // Correct mouse X coordinate calculation
        int x = static_cast<int>(mousePos.x / ((screenWidth - 200.0f) / gridSize));
        int y = static_cast<int>(mousePos.y / (screenHeight / gridSize));

        // Loop over a square region around the mouse position
        for (int offsetY = -brushRadius; offsetY <= brushRadius; ++offsetY)
        {
          for (int offsetX = -brushRadius; offsetX <= brushRadius; ++offsetX)
          {
            int neighborX = x + offsetX;
            int neighborY = y + offsetY;

            // Check if the neighboring point is within bounds
            if (neighborX >= 0 && neighborX < gridSize && neighborY >= 0 && neighborY < gridSize)
            {
              // Add density and velocity to the neighboring cell
              float distance = std::sqrt(offsetX * offsetX + offsetY * offsetY);
              float influence = std::max(0.0f, 1.0f - (distance / brushRadius));  // Influence decreases with distance

              fluid.add_density(neighborX, neighborY, densityAmount * influence);
              fluid.add_velocity(
                  neighborX, neighborY, GetMouseDelta().x * velocityFactor * influence, GetMouseDelta().y * velocityFactor * influence);
            }
          }
        }
      }
    }

    fluid.step();

    BeginDrawing();
    ClearBackground(BLACK);

    // Render fluid simulation
    fluid.render(screenWidth - 200, screenHeight, showVelocity);

    // Draw control panel
    DrawRectangle(screenWidth - 200, 0, 200, screenHeight, DARKGRAY);

    // Draw parameter selection and values
    for (int i = 0; i < paramCount; ++i)
    {
      Color textColor = (i == selectedParam) ? YELLOW : WHITE;
      DrawText(paramNames[i], screenWidth - 190, 20 + i * 40, 20, textColor);

      std::string valueStr;
      switch (i)
      {
        case 0:
          valueStr = float_to_string(diffusion);
          break;
        case 1:
          valueStr = float_to_string(viscosity);
          break;
        case 2:
          valueStr = float_to_string(timestep);
          break;
        case 3:
          valueStr = float_to_string(densityAmount);
          break;
        case 4:
          valueStr = float_to_string(velocityFactor);
          break;
        case 5:
          valueStr = std::to_string(brushRadius);
          break;
      }
      DrawText(valueStr.c_str(), screenWidth - 190, 45 + i * 40, 16, textColor);
    }

    // Draw instructions
    DrawText("Controls:", screenWidth - 190, screenHeight - 180, 20, WHITE);
    DrawText("Up/Down: Select Parameter", screenWidth - 190, screenHeight - 150, 16, LIGHTGRAY);
    DrawText("Left/Right: Adjust Value", screenWidth - 190, screenHeight - 130, 16, LIGHTGRAY);
    DrawText("V: Toggle Velocity/Density", screenWidth - 190, screenHeight - 110, 16, LIGHTGRAY);
    DrawText("Left Click: Add Fluid", screenWidth - 190, screenHeight - 90, 16, LIGHTGRAY);

    // Render current view mode
    std::string viewMode = showVelocity ? "Velocity View" : "Density View";
    DrawText(viewMode.c_str(), screenWidth - 190, screenHeight - 50, 20, WHITE);

    EndDrawing();
  }

  CloseWindow();
  return 0;
}
