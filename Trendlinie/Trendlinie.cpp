#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <stdexcept>
#include <cctype>
#include <cmath>
#include <algorithm>

// ================= CSV-Helfer =================
static inline std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b-1]))) --b;
    return s.substr(a, b - a);
}
static std::vector<std::string> split_semicolon(const std::string& s) {
    std::vector<std::string> cells; std::string cur;
    for (char ch : s) { if (ch == ';') { cells.push_back(trim(cur)); cur.clear(); } else { cur.push_back(ch); } }
    cells.push_back(trim(cur)); return cells;
}
static bool to_double(const std::string& in, double& out) {
    if (in.empty()) return false;
    std::string t; t.reserve(in.size());
    for (char ch : in) t.push_back(ch == ',' ? '.' : ch);
    try { size_t idx=0; out = std::stod(t, &idx); return idx == t.size(); } catch (...) { return false; }
}

// ================= Tabellen-Ausgabe =================
static void print_border(const std::vector<size_t>& w) {
    std::cout << '+'; for (size_t wi : w) { for (size_t i=0;i<wi+2;++i) std::cout << '-'; std::cout << '+'; } std::cout << '\n';
}
static void print_row(const std::vector<std::string>& row, const std::vector<size_t>& w) {
    std::cout << '|';
    for (size_t c=0;c<w.size();++c) {
        const std::string cell = c<row.size()? row[c] : "";
        std::cout << ' ' << cell; for (size_t i=cell.size(); i<w[c]; ++i) std::cout << ' '; std::cout << '|';
    }
    std::cout << '\n';
}

// ================= Regression =================
struct LinFit { double m=0, b=0, r2=0; size_t n=0; };
static LinFit linear_fit(const std::vector<double>& x, const std::vector<double>& y) {
    LinFit f; const size_t n=x.size(); if (n<2) return f;
    long double Sx=0,Sy=0,Sxx=0,Sxy=0,Syy=0;
    for (size_t i=0;i<n;++i){ Sx+=x[i]; Sy+=y[i]; Sxx+=x[i]*x[i]; Sxy+=x[i]*y[i]; Syy+=y[i]*y[i]; }
    const long double N=n, denom=N*Sxx - Sx*Sx; if (denom==0) return f;
    const long double m=(N*Sxy - Sx*Sy)/denom, b=(Sy - m*Sx)/N;
    const long double sst = Syy - (Sy*Sy)/N;
    long double sse=0; for (size_t i=0;i<n;++i){ const long double e=y[i] - (m*x[i]+b); sse+=e*e; }
    f.m=(double)m; f.b=(double)b; f.r2 = (sst>0)? (1.0L - sse/sst) : 0.0L; f.n=n; return f;
}

// ================= ASCII-Plot =================
struct Range { double lo, hi; };
static Range extend_range(double lo, double hi) {
    if (lo==hi) { const double eps = (std::abs(lo)+1.0)*0.05; lo-=eps; hi+=eps; }
    const double span = hi - lo; return { lo - 0.05*span, hi + 0.05*span };
}
static void plot_scatter_with_line(const std::vector<double>& X, const std::vector<double>& Y,
                                   const LinFit& F, const std::string& xlabel, const std::string& ylabel,
                                   int WIDTH=80, int HEIGHT=24)
{
    // Bereiche bestimmen (inkl. Regressionsgerade an Endpunkten)
    double xmin=*std::min_element(X.begin(),X.end());
    double xmax=*std::max_element(X.begin(),X.end());
    double ymin=*std::min_element(Y.begin(),Y.end());
    double ymax=*std::max_element(Y.begin(),Y.end());
    const double yL = F.m*xmin + F.b, yR = F.m*xmax + F.b;
    ymin = std::min(ymin, std::min(yL,yR));
    ymax = std::max(ymax, std::max(yL,yR));
    auto xr = extend_range(xmin, xmax);
    auto yr = extend_range(ymin, ymax);

    // Raster initialisieren
    std::vector<std::string> grid(HEIGHT, std::string(WIDTH, ' '));

    auto x_to_col = [&](double x)->int {
        double t = (x - xr.lo) / (xr.hi - xr.lo); t = std::clamp(t, 0.0, 1.0);
        return (int)std::llround(t * (WIDTH-1));
    };
    auto y_to_row = [&](double y)->int {
        double t = (y - yr.lo) / (yr.hi - yr.lo); t = std::clamp(t, 0.0, 1.0);
        return (int)std::llround((1.0 - t) * (HEIGHT-1)); // oben=0
    };

    auto put_char = [&](int r, int c, char ch) {
        if (r<0||r>=HEIGHT||c<0||c>=WIDTH) return;
        // Priorität: Datenpunkt 'x' > Linie '*' > Achsen '+', '-', '|'
        char& dst = grid[r][c];
        if (dst=='x') return;
        if (dst=='*' && ch!='x') return;
        dst = ch;
    };

    // Achsen (bei 0, falls im Bereich), sonst Rahmen
    bool drawY0 = (yr.lo<=0 && 0<=yr.hi);
    bool drawX0 = (xr.lo<=0 && 0<=xr.hi);
    int row0 = drawY0 ? y_to_row(0.0) : (HEIGHT-1);
    int col0 = drawX0 ? x_to_col(0.0) : 0;

    // X-Achse
    for (int c=0;c<WIDTH;++c) put_char(row0, c, '-');
    // Y-Achse
    for (int r=0;r<HEIGHT;++r) put_char(r, col0, '|');
    // Ursprung
    put_char(row0, col0, '+');

    // Regressionsgerade
    for (int c=0; c<WIDTH; ++c) {
        double x = xr.lo + (xr.hi - xr.lo) * (double)c / (double)(WIDTH-1);
        double y = F.m*x + F.b;
        int r = y_to_row(y);
        put_char(r, c, '*');
    }

    // Datenpunkte
    for (size_t i=0;i<X.size();++i) {
        int c = x_to_col(X[i]);
        int r = y_to_row(Y[i]);
        put_char(r, c, 'x');
    }

    // Ausgabe
    std::cout << "\nKoordinatensystem (x: links→rechts, y: unten→oben)\n";
    for (const auto& row : grid) std::cout << row << '\n';

    // Skalenhinweise
    std::cout << "x ~ [" << xr.lo << ", " << xr.hi << "]   "
              << "y ~ [" << yr.lo << ", " << yr.hi << "]\n";
    std::cout << "Legende: x = Datenpunkt, * = Regressionsgerade, +-| = Achsen\n";
}

// ================= Hauptprogramm =================
int main(int argc, char* argv[]) {
    if (argc < 2) { std::cerr << "Pfad zur CSV-Datei fehlt.\n"; return 1; }
    const char* path = argv[argc - 1];

    std::ifstream in(path);
    if (!in.is_open()) { std::cerr << "Fehler beim Öffnen: " << path << "\n"; return 1; }

    std::vector<std::vector<std::string>> rows;
    std::string line;
    while (std::getline(in, line)) rows.push_back(split_semicolon(line));
    in.close();
    if (rows.empty()) { std::cerr << "Datei ist leer.\n"; return 1; }

    // Spaltenzahl, Breiten
    size_t cols = 0; for (auto& r: rows) cols = std::max(cols, r.size());
    for (auto& r: rows) r.resize(cols);
    std::vector<size_t> w(cols,0);
    for (auto& r: rows) for (size_t c=0;c<cols;++c) w[c]=std::max(w[c], r[c].size());

    // Tabelle
    print_border(w); print_row(rows[0], w); print_border(w);
    for (size_t r=1;r<rows.size();++r) print_row(rows[r], w);
    print_border(w);

    // Spaltenübersicht
    std::cout << "\nSpalten (Index: Name):\n";
    for (size_t c=0;c<cols;++c)
        std::cout << "  [" << c << "] " << (rows[0][c].empty()? "(leer)" : rows[0][c]) << '\n';

    // Auswahl
    size_t ix, iy;
    std::cout << "\nIndex der X-Spalte: "; if (!(std::cin >> ix) || ix>=cols) { std::cerr << "Ungültig.\n"; return 1; }
    std::cout << "Index der Y-Spalte: "; if (!(std::cin >> iy) || iy>=cols) { std::cerr << "Ungültig.\n"; return 1; }

    // Daten extrahieren
    std::vector<double> X, Y; X.reserve(rows.size()); Y.reserve(rows.size());
    for (size_t r=1;r<rows.size();++r) {
        double x,y; if (to_double(rows[r][ix],x) && to_double(rows[r][iy],y)) { X.push_back(x); Y.push_back(y); }
    }
    if (X.size()<2) { std::cerr << "Zu wenige gültige Zahlenpaare.\n"; return 1; }

    // Fit
    LinFit F = linear_fit(X,Y);

    // Ausgabe Fit
    std::cout.setf(std::ios::fixed); std::cout.precision(6);
    std::cout << "\nTrendlinie (n=" << F.n << "): y = " << F.m << " * x + " << F.b << "\n";
    std::cout << "R^2 = " << F.r2 << "\n";

    // Plot
    plot_scatter_with_line(X, Y, F, rows[0][ix], rows[0][iy], 80, 24);
    return 0;
}
