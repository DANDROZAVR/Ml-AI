#include <bits/stdc++.h>
#include <Windows.h>

#define pb push_back
#define ld long double



ld sigmoid(ld x) {
  return 1 / (1 + exp(-x));
}
std::mt19937_64 gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());
class random {
public:
  ld segment11() const {
    ld x = (int) (gen() % 2001) - 1000;
    x /= 1000;
    return x;
  }
  inline ld segment01() const {
    return (ld) (gen() % 1001) / 1000;
  }
};

class Edge {
public:
  int from, to;
  ld val;
  Edge(int _from, int _to, ld _val) {
    from = _from;
    to = _to;
    val = _val;
  }
};


class NeuralNetwork {
  std::vector<int> neurons;
  std::vector<std::vector<int> > num;
  std::vector<ld> y, s, delta, WaSum;
  std::vector<std::vector<Edge> > edge;
  int N, sz;
  ld const1, const2;
public:
  void InitGen() {
    num.resize(neurons.size());
    edge.resize(neurons.size());
    WaSum.resize(neurons.size());
    num[0].resize(neurons[0]);
    random rand;
    for (int i = 0; i < neurons[0]; ++i)
      num[0][i] = N++;
    for (int lay = 0; lay + 1 < neurons.size(); ++lay) {
      num[lay + 1].resize(neurons[lay + 1]);
      for (int j = 0; j < neurons[lay + 1]; ++j)
        num[lay + 1][j] = N++;
      for (int i = 0; i < neurons[lay]; ++i)
        for (int j = 0; j < neurons[lay + 1]; ++j) {
          ld x = rand.segment11();
          edge[lay].pb(Edge(num[lay][i], num[lay + 1][j], x));
          ++sz;
        }
    }
    y.resize(N);
    s.resize(N);
    delta.resize(N);
    WaSum.resize(N);
  }
  void prec_gradient(std::vector<ld> &Test, ld res, int MinPrice, int MaxPrice) {
    fill(y.begin(), y.end(), 0);
    fill(s.begin(), s.end(), 0);
    for (int i = 0; i < neurons[0]; ++i) {
      s[num[0][i]] = 0;
      y[num[0][i]] = sigmoid(Test[i]);
    }
    for (int lay = 0; lay + 1 < neurons.size(); ++lay) {
      for (auto ed : edge[lay]) {
        s[ed.to] += y[ed.from] * ed.val;
      }
      for (int i = 0; i < neurons[lay + 1]; ++i) {
        y[num[lay + 1][i]] = sigmoid(s[num[lay + 1][i]]);
      }
    }
    const1 = res - MinPrice;
    const2 = MaxPrice - MinPrice;
  }
  void calc_gradient(std::vector<ld> &gradient) {
    fill(WaSum.begin(), WaSum.end(), 0);
    delta[N - 1] = -2 * const1 * const2 + const2 * const2 * y[N - 1] * 2;
    for (int lay = (int) (neurons.size()) - 2; lay >= 0; --lay) {
      for (auto i : edge[lay]) {
        WaSum[i.from] += i.val * delta[i.to];
      }
      for (int i = 0; i < neurons[lay]; ++i) {
        int v = num[lay][i];
        delta[v] = y[v] * (1 - y[v]) * WaSum[v];
      }
    }
    int id = 0;
    for (int lay = 0; lay + 1 < neurons.size(); ++lay)
      for (auto ed : edge[lay]) {
        gradient[id++] += delta[ed.to] * y[ed.from];
      }
  }
  void addNode(int _node) {
    neurons.pb(_node);
  }
  int getSize() const {
    return neurons.size();
  }
  void changeEdgesWeight(ld speed, std::vector < ld > & gradient) {
    int now = 0;
    for (int lay = 0; lay + 1 < neurons.size(); ++lay)
      for (int i = 0; i < edge[lay].size(); ++i) {
        edge[lay][i].val -= speed * gradient[now++];
      }
  }
  int getY(int index) const {
    return y[index];
  }
  ld getScore(ld &Last, int MinPrice, int MaxPrice, std::vector < std::vector < ld > > &test,
                        std::vector < int > &otv, int cnt_test, int Znam) {
    ld sum = 0;
    for (int i = 0; i < cnt_test; ++i) {
      prec_gradient(test[i], otv[i], MinPrice, MaxPrice);
      ld act = getY(N - 1);
      ld Price = act * (MaxPrice - MinPrice) + MinPrice;
      sum += (Price - otv[i]) * (Price - otv[i]);
    }
    Last = 1 - (sum / Znam);
    return 1 - (sum / Znam);
  }
};



const int MAXN = 3e5;
class solution {
  NeuralNetwork ml;
  std::vector< std::vector < ld > > test;
  std::vector<ld> eps, sigma;
  std::vector < int > otv;
  int cnt_test;
  double EpsAns, Znam;
  int n, MaxPrice = -1, MinPrice = 1e18, child;
  ld speed, koef1, beta, alpha;
public:
  void Init() {
    otv.resize(MAXN);
    test.resize(MAXN);
    cnt_test = 100000;
    n = 11;
    eps.resize(n);
    sigma.resize(n);
    for (int i = 0; i < cnt_test; ++i)
      test[i].resize(n);
  };
  void preGen() {
    Init();
    ReadData();
    NormalizationData();
    ml.addNode(n);
    ml.addNode(5);
    ml.addNode(5);
    ml.addNode(1);
    ml.InitGen();
  }
  void run() {
    ld Last, PrLast = -1e9;
    std::vector<ld> gradient(ml.getSize());
    int Iter = 0;
    speed = 0.01;
    int clc = 1;
    while (true) {
      std::cerr << "SCORE = " << ++Iter << " " << ml.getScore(Last, MinPrice, MaxPrice, test, otv, cnt_test, Znam) << " " << speed << std::endl;
      fill(gradient.begin(), gradient.end(), 0);
      int cnt = 0;
      for (int i = 0; i < cnt_test; ++i) {
        ml.prec_gradient(test[i], otv[i], MinPrice, MaxPrice);
        ml.calc_gradient(gradient);
        ++cnt;
      }
      for (int i = 0; i < ml.getSize(); ++i)
        gradient[i] /= cnt;
      ld mxx = -1e18;
      for (int i = 0; i < ml.getSize(); ++i) {
        mxx = std::max(mxx, gradient[i]);
      }
      if (abs(mxx) < 0.000000000001) break;
      for (int i = 0; i < ml.getSize(); ++i) {
        gradient[i] /= abs(mxx);
      }
      int Now = 0;
      ml.changeEdgesWeight(speed, gradient);
      if (Last > PrLast * koef1)
        speed *= beta;
      else
        speed *= alpha;
      PrLast = Last;
    }
  }
  ld Stod(std::string s){
    ld ans = 0;
    while(s.size() && s[0] != '.'){
      ans = ans * 10 + s[0] - '0';
      s.erase(0, 1);
    }
    if (s.size()){
      while(s.size() > 1){
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
    std::ifstream in("NewTest2.csv"); // local input file from
    std::string s;
    getline(in, s);
    std::vector<std::string> v;
    for (int i = 0; i < cnt_test; ++i) {
      getline(in, s);
      v.pb(s);
    }
//    shuffle(v.begin(), v.end(), gen);
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
          MaxPrice = std::max(MaxPrice, otv[num]);
          MinPrice = std::min(MinPrice, otv[num]);
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
    memset(Eps, 0, sizeof(Eps) * n);
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

  void NormalizationData() { // decrease
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
};




int main() {
  setlocale(LC_ALL, "Russian");
  SetConsoleCP(1251);
  SetConsoleOutputCP(1251);

  solution main;
  main.Init();
  main.ReadData();
  main.preGen();
  main.NormalizationData();
  main.run();
  return 0;
}


