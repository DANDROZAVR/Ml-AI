#include <bits/stdc++.h>
#include <Windows.h>
#define pb push_back
#define ld long double
#define S second
#define F first
#define vvld std::vector < std::vector < ld > >

std::mt19937_64 gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());
ld sigmoid(ld x) {
  return 1 / (1 + exp(x));
}
class random {
public:
  ld segment11() {
    ld x = (int) (gen() % 2001) - 1000;
    x /= 1000;
    return x;
  }

  inline ld segment01() {
    return (ld) (gen() % 1001) / 1000;
  }
} random;

class genes {
  ld BEST = 0.0646384;
  std::vector<ld> OTV;
  ld b1 = 5;
  ld b2 = 0;
  ld b3 = -3;
  ld b4 = 6.8;
public:
  ld Cost(std::vector<ld> &genom, int cnt_test, int n, std::vector < std::vector < ld > > &test, std::vector < int > otv, int MinPrice, int std::maxPrice) {
    ld ans = 0;
    for (int i = 0; i < cnt_test; ++i) {
      ld sum1 = 0, sum2 = 0, sum3 = 0;
      for (int j = 0; j < n; ++j)
        sum1 += genom[j] * test[i][j];
      for (int j = n; j < n * 2; ++j)
        sum2 += genom[j] * test[i][j - n];
      for (int j = n * 2; j < n * 3; ++j)
        sum3 += genom[j] * test[i][j - n * 2];
      ld activ1 = sigmoid(sum1 + b1);
      ld activ2 = sigmoid(sum2 + b2);
      ld activ3 = sigmoid(sum3 + b3);

      ld Final = activ1 * genom[n * 3];
      Final += activ2 * genom[n * 3 + 1];
      Final += activ3 * genom[n * 3 + 2];
      ld Factive = sigmoid(Final + b4);
      ld Price = Factive * (std::maxPrice - MinPrice) + MinPrice;
      ans += (Price - otv[i]) * (Price - otv[i]);
    }
    return ans;
  }

  ld Score(std::vector<ld> &genom, int Znam, int cnt_test, int n, std::vector < std::vector < ld > > &test, std::vector < int > otv, int MinPrice, int std::maxPrice) {
    ld ans = 1 - (Cost(genom, cnt_test, n, test, otv, MinPrice, std::maxPrice) / Znam);
    if (ans > BEST) {
      BEST = ans;
      ofstream out("output.txt");
      out << ans << " ";
      for (auto j : genom) out << j << " ";
      out << endl;
      out.close();
    }
    return ans;
  }
  vvld evolution(vvld genom, int N, int Znam, int PopSz, int child, int cnt_test, int n, std::vector < std::vector < ld > > &test, std::vector < int > otv, int MinPrice, int std::maxPrice) {
    vvld ne;
    int m = genom.size();
    std::pair < ld, int > cost[m];
    for (int i = 0; i < m; ++i){
      cost[i].F = Cost(genom[i], cnt_test, n, test, otv, MinPrice, std::maxPrice);
      cost[i].S = i;
    }
    sort(cost, cost + m, [&](std::pair < ld, int > a, std::pair < ld, int > b){
      return a.F < b.F;
    });
    //std::cer<<fixed<<setprecision(7)<<Score(genom[cost[0].S], Znam, cnt_test, n, test, otv, MinPrice, std::maxPrice)<<endl; // progress genes
    m = PopSz;
    while(ne.size() < PopSz){
      for (int kol = 0; kol < child; ++kol)
        for (int j = 0; j < m; ++j){
          auto NewGen = genom[cost[j].S];
          if (gen() % 4){
            for (int i = 0; i < N; ++i)
              if (gen() % 27 == 0){
                NewGen[i] = random.segment11();
              }else
              if (gen() % 7 == 0){
                NewGen[i] += random.segment01() / 10;
                NewGen[i] = min(NewGen[i], (ld)(1));
                NewGen[i] = std::max(NewGen[i], (ld)(-1));
              }
          }
          ne.pb(NewGen);
        }
    }
    return ne;
  }
};

class solution {
  int PopSz, cnt_test, n, std::maxPrice = -1, MinPrice = 1e18, child, n2, N;
  double Mut, Pnxt, ch_by_iter, EpsAns, Znam;
  const int std::maxN = 3e5;
  std::vector < std::vector<ld> > test;
  std::vector<ld> eps, sigma;
  std::vector < int > otv;
public:
  int getN() const {
    return N;
  }
  void Init(int _PopSz = 140, double _Mut = 0.5, int _child = 1, int _cnt_test = 140000, int _n = 11, int _n2 = 3, ld _ch_by_iter = 0.6) {
    test.resize(std::maxN);
    otv.resize(std::maxN);
    PopSz = _PopSz;
    Mut = _Mut;
    child = _child;
    cnt_test = _cnt_test;
    n = _n;
    n2 = _n2;
    N = n2 + n * n2;
    ch_by_iter = _ch_by_iter;
    eps.resize(n);
    sigma.resize(n);
    for (int i = 0; i < cnt_test; ++i)
      test[i].resize(n);
  }
  ld Stod(std::string s) {
    ld ans = 0;
    while (s.size() && s[0] != '.') {
      ans = ans * 10 + s[0] - '0';
      s.erase(0, 1);
    }
    if (s.size()) {
      while (s.size() > 1) {
        ld x = s.back() - '0';
        ans += x;
        ans /= 10;
        s.pop_back();
      }
    }
    return ans;
  }
  void ReadData() {
    std::string z = "5";
    ifstream in("NewTest.csv");
    std::string s;
    getline(in, s);
    std::vector<std::string> v;
    for (int i = 0; i < cnt_test; ++i) {
      getline(in, s);
      v.pb(s);
    }
    shuffle(v.begin(), v.end(), gen);

    int num = 0;
    for (std::string &s : v) {
      int i = 0;
      int cnt = 0, now = 0;
      while (i < s.size()) {
        int Pi = i;
        while (i < s.size() && s[i] != ',') ++i;
        std::string q = s.substr(Pi, i - Pi);
        if (cnt == 11) {
          otv[num] = stod(q);
          EpsAns += otv[num];
          std::maxPrice = std::max(std::maxPrice, otv[num]);
          MinPrice = min(MinPrice, otv[num]);
        } else if ((cnt >= 3 && cnt <= 5) || (cnt >= 8 && cnt <= 10)) {
          if (q == "nan" || q.empty()) {
            test[num][now] = -1;
          } else
            test[num][now] = Stod(q);
        } else if (cnt == 6) {
          if (q == "дизель") {
            test[num][now] = 1;
            ++now;
            test[num][now] = 0;
            ++now;
            test[num][now] = 0;
          } else if (q == "бензин") {
            test[num][now] = 0;
            ++now;
            test[num][now] = 1;
            ++now;
            test[num][now] = 0;
          } else {
            test[num][now] = 0;
            ++now;
            test[num][now] = 0;
            ++now;
            test[num][now] = 1;
          }
        } else if (cnt == 7) {
          if (q == "авто") {
            test[num][now] = 1;
            ++now;
            test[num][now] = 0;
          } else {
            test[num][now] = 0;
            ++now;
            test[num][now] = 1;
          }
        } else --now;
        ++cnt;
        ++now;
        ++i;
      }
      ++num;
    }
    int Eps[n];
    fill(Eps, Eps + n, 0);
    for (int i = 0; i < cnt_test; ++i) {
      for (int j = 0; j < n; ++j)
        if (test[i][j] != -1) {
          Eps[j] += test[i][j];
        }
    }
    for (int i = 0; i < cnt_test; ++i) {
      for (int j = 0; j < n; ++j)
        if (test[i][j] == -1) {
          test[i][j] = Eps[j];
        }
    }

    EpsAns /= cnt_test;
    for (int i = 0; i < cnt_test; ++i)
      Znam += (EpsAns - otv[i]) * (EpsAns - otv[i]);
  }
  void NormData() {
    for (int i = 0; i < cnt_test; ++i)
      for (int j = 0; j < n; ++j)
        eps[j] += test[i][j];
    for (int j = 0; j < n; ++j) {
      eps[j] /= cnt_test;
      for (int i = 0; i < cnt_test; ++i)
        sigma[j] += (test[i][j] - eps[j]) * (test[i][j] - eps[j]);
      sigma[j] /= cnt_test - 1;
      for (int i = 0; i < cnt_test; ++i)
        test[i][j] = (test[i][j] - eps[j]) / sigma[j];
    }
  }
  void run() {
    genes population;
    std::vector < std::vector < ld > > allgenes;
    for (int j = 0; j < 1000; ++j){
      std::vector < ld > v;
      for (int i = 0; i < N; ++i)
        v.pb(random.segment11());
      allgenes.pb(v);
    }

    for (int iter = 0; iter < 100000; ++iter){
      auto V = population.evolution(allgenes, N, Znam, PopSz, child, cnt_test, n, test, otv, MinPrice, std::maxPrice);
      allgenes = V;
    }
    // ready for counting the data, taking into account the genes
  }
};




int main(){
  setlocale(LC_ALL,"Russian");
  SetConsoleCP(1251);
  SetConsoleOutputCP(1251);

  solution solut;
  solut.Init();
  solut.ReadData();
  solut.NormData();
  solut.run();
  return 0;
}

