#include <iostream>
#include <cstring>
#include <ncurses.h>
#include <vector>       // For std::vector
#include <thread>       // For std::thread, std::thread::hardware_concurrency
#include <cmath>        // For std::ceil, std::min
#include <cstddef>      // For size_t
#include <stdexcept>    // For std::runtime_error (optional error checking)
			//
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

// --- Multithreaded Mandelbrot Function ---
template <typename float_t>
void mandelbrot_threaded(MandelbrotData *mandel_data, float_t xmid, float_t xhalfrange, float_t ymid, float_t yhalfrange, unsigned int num_threads) {
    if (num_threads == 0) {
        throw std::runtime_error("Number of threads must be greater than zero.");
    }
     if (mandel_data->x_res <= 0 || mandel_data->y_res <= 0) {
        throw std::runtime_error("Invalid dimensions in MandelbrotData.");
    }


    const float_t x1 = xmid - xhalfrange;
    const float_t x2 = xmid + xhalfrange;
    const float_t y1 = ymid - yhalfrange;
    const float_t y2 = ymid + yhalfrange;

    const int x_res = mandel_data->x_res;
    const int y_res = mandel_data->y_res;

    // Calculate the step size between pixels in the complex plane
    // Avoid division by zero if resolution is 1
    const float_t xRes0 = (x_res > 1) ? (x2 - x1) / (x_res - 1) : 0;
    const float_t yRes0 = (y_res > 1) ? (y2 - y1) / (y_res - 1) : 0;

    // --- Threading Logic ---
    std::vector<std::thread> threads;
    threads.reserve(num_threads); // Pre-allocate space for threads

    // Calculate the approximate number of rows each thread should process
    // Use ceil to ensure all rows are covered
    const size_t total_rows = static_cast<size_t>(y_res);
    const size_t rows_per_thread = static_cast<size_t>(std::ceil(static_cast<double>(total_rows) / num_threads));

    // --- Launch Threads ---
    // Each thread calculates a horizontal strip of the Mandelbrot set
    for (unsigned int i = 0; i < num_threads; ++i) {
        // Determine the start and end row for this thread
        const size_t py_start = std::min(i * rows_per_thread, total_rows);
        const size_t py_end = std::min(py_start + rows_per_thread, total_rows);

        // Skip creating a thread if the range is empty (can happen if num_threads > y_res)
        if (py_start >= py_end) {
            continue;
        }

        // Create and launch a thread using a lambda function
        threads.emplace_back([=]() { // Capture necessary variables by value
            // Per-thread variables for Mandelbrot calculation
            float_t x0, y0, x_, y_, x, y, sq_x, sq_y;

            // Iterate over the rows assigned to this thread
            for (size_t py = py_start; py < py_end; ++py) {
                 // Calculate the initial imaginary part (y0) for this row
                 y0 = py * yRes0;
                 y0 += y1;

                 // Iterate over all columns (px) for the current row
                for (size_t px = 0; px < static_cast<size_t>(x_res); ++px) {
                    // Reset variables for each pixel
                    x = y = 0;
                    sq_x = 0;
                    sq_y = 0;

                    // Calculate the initial real part (x0) for this pixel
                    x0 = px * xRes0;
                    x0 += x1;

                    int n = 0; // Iteration counter

                    // The core Mandelbrot iteration loop
                    // Check if the point escapes the circle |z| < 2 (i.e., sq_x + sq_y < 4)
                    // or if the maximum number of iterations is reached.
                    for (n = 0; ((sq_x + sq_y) < 4.0) && (n < MaxIter); ++n) {
                        // Calculate the next iteration: z = z^2 + c
                        // where z = x + yi, c = x0 + y0i
                        // z^2 = (x+yi)^2 = x^2 - y^2 + 2xyi
                        // Next x: (x^2 - y^2) + x0
                        // Next y: (2xy) + y0
                        y_ = (2.0 * x * y) + y0; // Calculate new y first using old x and y
                        x_ = (sq_x - sq_y) + x0; // Calculate new x using old sq_x and sq_y

                        // Update x and y for the next iteration
                        x = x_;
                        y = y_;

                        // Pre-calculate squares for the escape condition and next iteration
                        sq_x = x * x;
                        sq_y = y * y;
                    }

                    // Store the result (number of iterations) in the MandelbrotData structure
                    // This assumes mandel_data->set_pixel_data is thread-safe for different (px, py)
                    mandel_data->set_pixel_data(px, py, n);
                } // End loop over columns (px)
            } // End loop over rows (py)
        }); // End of lambda function and thread creation
    } // End loop for launching threads

    // --- Join Threads ---
    // Wait for all launched threads to finish their execution
    for (auto& t : threads) {
        if (t.joinable()) { // Check if the thread is joinable (i.e., not already joined)
            t.join();
        }
    }
}

int main(int argc, char **argv) {
    // Initialize ncurses
    initscr();
    keypad(stdscr, TRUE);
    noecho();
    cbreak();

    // Get the number of hardware threads available (returns 0 if undetectable)
    unsigned int hardware_threads = std::thread::hardware_concurrency();

    // Use the hardware thread count, or default to 2 if detection fails
    unsigned int num_threads_to_use = (hardware_threads == 0) ? 2 : hardware_threads;

    printw("ASCII Mandel 1.0\n\n");
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
      mandelbrot_threaded<double>(&mandel_data, xmid, xhalfrange, ymid, xhalfrange / aspectratio, num_threads_to_use);

      mandel_data.normalize();
      mandel_data.print();

      printw("X: %f Y: %f zoom: %fx N: %d threads: %d", xmid, ymid, 2.5 / xhalfrange, MaxIter, num_threads_to_use);
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
