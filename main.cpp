#include <iostream>
#include <cstring>
#include <ncurses.h>

size_t MaxIter = 500;
const char* Pallete = " .'`^\",:;Il!i><~+_-?][}{1)(|\\/tfjrxnuvczXYUJCLQ0OZmwqpdbkhao*#MW&8%B@$";

// MandelbrotData stores pixel data for each frame as a floating-point
// in the range of [0, 1] based on iterations spent on each pixel
struct MandelbrotData {
    int x_res, y_res;
    double* data;

    MandelbrotData() = delete;

    // Default constructor, allocate memory
    MandelbrotData(int xres, int yres) : x_res(xres), y_res(yres)
    {
        data = new double[x_res * y_res];
        memset(data, 0, size_t(x_res) * size_t(y_res) * sizeof(double));
    }

    ~MandelbrotData()
    {
        delete[] data;
    }

    // Set pixel data based on no. of iterations
    void set_pixel_data(int x, int y, size_t n)
    {
        *(data + y*x_res + x) = double(n) / double(MaxIter);
    }

    double get_pixel_data(int x, int y)
    {
        return *(data + y*x_res + x);
    }

    // Map all of the pixel values between [0.0, 1.0]
    void normalize()
    {
      size_t i;
      double minval = 1, maxval = 0;

      for (i = 0; i < x_res * y_res; ++i) {
        if (data[i] > maxval)
          maxval = data[i];
        if (data[i] < minval)
          minval = data[i];
      }

      if (maxval == minval)
        return;

      double norm_factor = 1. / (maxval - minval);

      for (i = 0; i < x_res * y_res; ++i)
        data[i] = (data[i] - minval) * norm_factor;
    }

    void print()
    {
        move(0, 0);

        const int PalleteLen = strlen(Pallete);
        for (size_t py = 0; py < y_res; py++)
        {
            for (size_t px = 0; px < x_res; px++)
                mvaddch(py, px, Pallete[int(get_pixel_data(px, py) * PalleteLen)]);
        }
        refresh();
    }
};

template <typename float_t>
void mandelbrot(MandelbrotData *mandel_data, float_t xmid, float_t xhalfrange, float_t ymid, float_t yhalfrange) {
  float_t x1, x2, y1, y2;

  int x_res, y_res;
  x_res = mandel_data->x_res; y_res = mandel_data->y_res;

  x1 = xmid - xhalfrange;
  x2 = xmid + xhalfrange;
  y1 = ymid - yhalfrange;
  y2 = ymid + yhalfrange;

  float_t xRes0 = (x2 - x1) / (x_res - 1);
  float_t yRes0 = (y2 - y1) / (y_res - 1);

  float_t x0, y0, x_, y_, x, y, sq_x, sq_y;

  for (size_t px = 0; px < x_res; ++px) {
    for (size_t py = 0; py < y_res; ++py) {
      x = y = 0;

      x0 = px * xRes0;
      x0 += x1;
      y0 = py * yRes0;
      y0 += y1;

      int n = 0;
      sq_x = x * x;
      sq_y = y * y;

      for (n = 0; ((sq_x + sq_y) < 4) && (n < MaxIter); ++n) {
        x_ = (sq_x - sq_y) + x0;
        y_ = (2 * x * y) + y0;

        x = x_;
        y = y_;
        sq_x = x * x;
        sq_y = y * y;
      }

      mandel_data->set_pixel_data(px, py, n);
    }
  }
}

int main(int argc, char **argv) {
    // Initialize ncurses
    initscr();
    keypad(stdscr, TRUE);
    noecho();
    cbreak();

    printw("ASCII Mandel 0.0\n\n");
    refresh();

    int x_res, y_res;
    getmaxyx(stdscr, y_res, x_res);

    attron(A_BOLD);
    printw("w/a/s/d");
    attroff(A_BOLD);
    printw(" to navigate\n");

    attron(A_BOLD);
    printw("Arrow keys");
    attroff(A_BOLD);
    printw(" to zoom in/out\n");

    attron(A_BOLD);
    printw("[ ]");
    attroff(A_BOLD);
    printw(" to decrease/increase contrast\n");

    attron(A_BOLD);
    printw("q");
    attroff(A_BOLD);
    printw(" to quit\n\n");

    refresh();

    printw("Press any key to continue...");
    refresh();
    getch();

    MandelbrotData mandel_data(x_res, y_res - 1); // Save one line for info

    double xmid = -0.75, ymid = 0;
    double xhalfrange = 2.5;
    double aspectratio = 1.5;
    while (true) {
      clear();
      mandelbrot<double>(&mandel_data, xmid, xhalfrange, ymid, xhalfrange / aspectratio);

      mandel_data.normalize();
      mandel_data.print();

      printw("X: %f Y: %f zoom: %fx N: %d", xmid, ymid, 2.5 / xhalfrange, MaxIter);
      refresh();
      int ch = getch();

      if (std::tolower(ch) == 'w') // Navigate up
        ymid -= 2 * (xhalfrange / aspectratio) / y_res;
      else if (std::tolower(ch) == 's') // Navigate down
        ymid += 2 * (xhalfrange / aspectratio) / y_res;
      else if (std::tolower(ch) == 'a') // Navigate left
        xmid -= 2 * xhalfrange / x_res;
      else if (std::tolower(ch) == 'd') // Navigate right
        xmid += 2 * xhalfrange / x_res;
      else if (ch == KEY_UP || ch == '+') // ZOOM IN
        xhalfrange -= 2 * xhalfrange / x_res;
      else if (ch == KEY_DOWN || ch == '-') // ZOOM OUT
        xhalfrange += 2 * xhalfrange / x_res;
      else if (std::tolower(ch) == ']') // Increase number of iterations
        MaxIter *= 1.1;
      else if (std::tolower(ch) == '[') // Decrease number of iterations
        MaxIter /= 1.1;
      else if (std::tolower(ch) == 'q') // Exit
        break;
    }

    endwin();
    return 0;
}
