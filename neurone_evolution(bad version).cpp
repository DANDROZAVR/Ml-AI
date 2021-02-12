#include <bits/stdc++.h>
#include <Windows.h>
#define pb push_back
#define ld long double
#define S second
#define F first
#define vvld vector < vector < ld > >

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
using namespace std;
class genes {
  ld BEST = 0.0646384;
  vector<ld> OTV;
  ld b1 = 5;
  ld b2 = 0;
  ld b3 = -3;
  ld b4 = 6.8;
public:
  ld Cost(vector<ld> &genom, int cnt_test, int n, vector < vector < ld > > &test, vector < int > otv, int MinPrice, int MaxPrice) {
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
      ld Price = Factive * (MaxPrice - MinPrice) + MinPrice;
      ans += (Price - otv[i]) * (Price - otv[i]);
    }
    return ans;
  }

  ld Score(vector<ld> &genom, int Znam, int cnt_test, int n, vector < vector < ld > > &test, vector < int > otv, int MinPrice, int MaxPrice) {
    ld ans = 1 - (Cost(genom, cnt_test, n, test, otv, MinPrice, MaxPrice) / Znam);
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
  vvld evolution(vvld genom, int N, int Znam, int PopSz, int child, int cnt_test, int n, vector < vector < ld > > &test, vector < int > otv, int MinPrice, int MaxPrice) {
    vvld ne;
    int m = genom.size();
    pair < ld, int > cost[m];
    for (int i = 0; i < m; ++i){
      cost[i].F = Cost(genom[i], cnt_test, n, test, otv, MinPrice, MaxPrice);
      cost[i].S = i;
    }
    sort(cost, cost + m, [&](pair < ld, int > a, pair < ld, int > b){
      return a.F < b.F;
    });
    //cerr<<fixed<<setprecision(7)<<Score(genom[cost[0].S], Znam, cnt_test, n, test, otv, MinPrice, MaxPrice)<<endl; // progress genes
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
                NewGen[i] = max(NewGen[i], (ld)(-1));
              }
          }
          ne.pb(NewGen);
        }
    }
    return ne;
  }
};

class solution {
  int PopSz, cnt_test, n, MaxPrice = -1, MinPrice = 1e18, child, n2, N;
  double Mut, Pnxt, ch_by_iter, EpsAns, Znam;
  const int MAXN = 3e5;
  vector < vector<ld> > test;
  vector<ld> eps, sigma;
  vector < int > otv;
public:
  int getN() const {
    return N;
  }
  void Init(int _PopSz = 140, double _Mut = 0.5, int _child = 1, int _cnt_test = 140000, int _n = 11, int _n2 = 3, ld _ch_by_iter = 0.6) {
    test.resize(MAXN);
    otv.resize(MAXN);
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
  ld Stod(string s) {
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
    string z = "5";
    ifstream in("NewTest.csv");
    string s;
    getline(in, s);
    vector<string> v;
    for (int i = 0; i < cnt_test; ++i) {
      getline(in, s);
      v.pb(s);
    }
    shuffle(v.begin(), v.end(), gen);

    int num = 0;
    for (string &s : v) {
      int i = 0;
      int cnt = 0, now = 0;
      while (i < s.size()) {
        int Pi = i;
        while (i < s.size() && s[i] != ',') ++i;
        string q = s.substr(Pi, i - Pi);
        if (cnt == 11) {
          otv[num] = stod(q);
          EpsAns += otv[num];
          MaxPrice = max(MaxPrice, otv[num]);
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
    vector < vector < ld > > allgenes;
    for (int j = 0; j < 1000; ++j){
      vector < ld > v;
      for (int i = 0; i < N; ++i)
        v.pb(random.segment11());
      allgenes.pb(v);
    }

    for (int iter = 0; iter < 100000; ++iter){
      auto V = population.evolution(allgenes, N, Znam, PopSz, child, cnt_test, n, test, otv, MinPrice, MaxPrice);
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

