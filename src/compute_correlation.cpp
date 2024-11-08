#include <bits/stdc++.h>
using namespace std;

struct WigEntry {
  string chrom;
  int pos;
  double value;
};

// Parse a single WIG fixed step header line
pair<string, int> parseHeader(const string &line) {
  string chrom = "";
  int start = -1;

  size_t chrom_pos = line.find("chrom=");
  size_t start_pos = line.find("start=");

  if (chrom_pos != string::npos) {
    size_t chrom_end = line.find(" ", chrom_pos);
    if (chrom_end == string::npos)
      chrom_end = line.length();
    chrom = line.substr(chrom_pos + 6, chrom_end - (chrom_pos + 6));
  }

  if (start_pos != string::npos) {
    size_t start_end = line.find(" ", start_pos);
    if (start_end == string::npos)
      start_end = line.length();
    string start_str = line.substr(start_pos + 6, start_end - (start_pos + 6));
    start = stoi(start_str);
  }

  return {chrom, start};
}

vector<WigEntry> readWig(const string &filename) {
  ifstream fin(filename);
  if (!fin.is_open()) {
    cerr << "Cannot open file: " << filename << endl;
    exit(1);
  }

  vector<WigEntry> entries;
  string line, curr_chrom;
  int curr_pos = 0;

  while (getline(fin, line)) {
    if (line.empty())
      continue;

    if (line.substr(0, 9) == "fixedStep") {
      auto [chrom, start] = parseHeader(line);
      if (chrom.empty() || start == -1) {
        cerr << "Invalid header line: " << line << endl;
        exit(1);
      }
      curr_chrom = chrom;
      curr_pos = start;
      continue;
    }

    try {
      double value = stod(line);
      entries.push_back({curr_chrom, curr_pos++, value});
      // cout << "curr_chrom" << curr_chrom << endl;
    } catch (const std::exception &e) {
      cerr << "Error parsing value: " << line << endl;
      exit(1);
    }
  }

  return entries;
}

double computeCorrelation(const vector<WigEntry> &wig1,
                          const vector<WigEntry> &wig2) {
  map<pair<string, int>, double> vals1, vals2;
  for (const auto &entry : wig1) {
    vals1[{entry.chrom, entry.pos}] = entry.value;
  }
  for (const auto &entry : wig2) {
    vals2[{entry.chrom, entry.pos}] = entry.value;
  }

  vector<double> x, y;
  for (const auto &[key, val1] : vals1) {
    auto it = vals2.find(key);
    if (it != vals2.end()) {
      x.push_back(val1);
      y.push_back(it->second);
    }
  }

  int n = x.size();
  if (n == 0) {
    cerr << "No overlapping positions found" << endl;
    return 0.0;
  }

  double sum_x = 0, sum_y = 0;
  for (int i = 0; i < n; i++) {
    sum_x += x[i];
    sum_y += y[i];
  }
  double mean_x = sum_x / n;
  double mean_y = sum_y / n;

  double numerator = 0;
  double denom_x = 0, denom_y = 0;

  for (int i = 0; i < n; i++) {
    double dx = x[i] - mean_x;
    double dy = y[i] - mean_y;
    numerator += dx * dy;
    denom_x += dx * dx;
    denom_y += dy * dy;
  }

  if (denom_x == 0 || denom_y == 0)
    return 0.0;
  return numerator / sqrt(denom_x * denom_y);
}

int main(int argc, char *argv[]) {
  ios_base::sync_with_stdio(false);
  cin.tie(nullptr);

  if (argc != 3) {
    cerr << "Usage: " << argv[0] << " <wig1> <wig2>" << endl;
    return 1;
  }

  // Add some basic logging
  cerr << "Reading file 1: " << argv[1] << endl;
  auto wig1 = readWig(argv[1]);
  cerr << "Reading file 2: " << argv[2] << endl;
  auto wig2 = readWig(argv[2]);

  cerr << "Computing correlation..." << endl;
  cerr << "File 1 entries: " << wig1.size() << endl;
  cerr << "File 2 entries: " << wig2.size() << endl;

  double correlation = computeCorrelation(wig1, wig2);
  cout << fixed << setprecision(3) << correlation << endl;

  return 0;
}
