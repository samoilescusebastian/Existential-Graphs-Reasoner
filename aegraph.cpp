// Copyright 2019 Luca Istrate, Danut Matei
#ifndef _AEGRAPH_H_
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <utility>
#include <cassert>
#include "./aegraph.h"
std::string strip(std::string s) {
  // deletes whitespace from the beginning and end of the string
  s.erase(0, s.find_first_not_of(" \n\r\t"));
  s.erase(s.find_last_not_of(" \n\r\t") + 1);
  return s;
}

std::pair<std::string, std::string> split_first(std::string s,
                                                char delimiter = ',') {
  // returns a pair that contains: <first_cut, rest_of_graph>

  int numOpen = 0;

  int len = s.size();
  for (int i = 0; i < len; i++) {
    char c = s[i];
    if (c == delimiter && numOpen == 0)
      return std::make_pair(strip(s.substr(0, i)), strip(s.substr(i + 1, len)));
    if (c == '[')
      numOpen += 1;
    if (c == ']')
      numOpen -= 1;
  }

  return {strip(s), std::string()};
}

std::vector<std::string> split_level(std::string s, char delimiter = ',') {
  // splits 's' into separate entities (atoms, subgraphs)

  std::vector<std::string> result;
  auto snd = s;
  while (true) {
    auto p = split_first(snd, delimiter);
    auto fst = p.first;
    snd = p.second;

    result.push_back(fst);
    if (snd.empty())
      return result;
  }
}

int AEGraph::num_subgraphs() const {
  return subgraphs.size();
}

int AEGraph::num_atoms() const {
  return atoms.size();
}

int AEGraph::size() const {
  return num_atoms() + num_subgraphs();
}

bool AEGraph::operator<(const AEGraph& other) const {
  return this->repr() < other.repr();
}

bool AEGraph::operator==(const AEGraph& other) const {
  return this->repr() == other.repr();
}

bool AEGraph::operator!=(const AEGraph& other) const {
  return this->repr() != other.repr();
}

AEGraph AEGraph::operator[](const int index) const {
  // offers an easier way of accessing the nested graphs
  if (index < num_subgraphs()) {
    return subgraphs[index];
  }

  if (index < num_subgraphs() + num_atoms()) {
    return AEGraph('(' + atoms[index - num_subgraphs()] + ')');
  }

  return AEGraph("()");
}

std::ostream& operator<<(std::ostream& out, const AEGraph& g) {
  out << g.repr();
  return out;
}

AEGraph::AEGraph(std::string representation) {
  // constructor that creates an AEGraph structure from a
  // serialized representation
  char left_sep = representation[0];
  char right_sep = representation[representation.size() - 1];

  assert((left_sep == '(' && right_sep == ')') ||
         (left_sep == '[' && right_sep == ']'));

  // if the left separator is '(' then the AEGraph is the entire
  // sheet of assertion
  if (left_sep == '(') {
    is_SA = true;
  } else {
    is_SA = false;
  }

  // eliminate the first pair of [] or ()
  representation = representation.substr(1, representation.size() - 2);

  // split the graph into separate elements
  auto v = split_level(representation);
  // add them to the corresponding vector
  for (auto s : v) {
    if (s[0] != '[') {
      atoms.push_back(s);
    } else {
      subgraphs.push_back(AEGraph(s));
    }
  }

  // also internally sort the new graph
  this->sort();
}

std::string AEGraph::repr() const {
  // returns the serialized representation of the AEGraph

  std::string left, right;
  if (is_SA) {
    left = '(';
    right = ')';
  } else {
    left = '[';
    right = ']';
  }

  std::string result;
  for (auto subgraph : subgraphs) {
    result += subgraph.repr() + ", ";
  }

  int len = atoms.size();
  if (len != 0) {
    for (int i = 0; i < len - 1; i++) {
      result += atoms[i] + ", ";
    }
    result += atoms[len - 1];
  } else {
    if (subgraphs.size() != 0)
      return left + result.substr(0, result.size() - 2) + right;
  }

  return left + result + right;
}

void AEGraph::sort() {
  std::sort(atoms.begin(), atoms.end());

  for (auto& sg : subgraphs) {
    sg.sort();
  }

  std::sort(subgraphs.begin(), subgraphs.end());
}

bool AEGraph::contains(const std::string other) const {
  // checks if an atom is in a graph
  if (find(atoms.begin(), atoms.end(), other) != atoms.end())
    return true;

  for (const auto& sg : subgraphs)
    if (sg.contains(other))
      return true;

  return false;
}

bool AEGraph::contains(const AEGraph& other) const {
  // checks if a subgraph is in a graph
  if (find(subgraphs.begin(), subgraphs.end(), other) != subgraphs.end())
    return true;

  for (const auto& sg : subgraphs)
    if (sg.contains(other))
      return true;

  return false;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(
    const std::string other) const {
  // returns all paths in the tree that lead to an atom like <other>
  std::vector<std::vector<int>> paths;

  int len_atoms = num_atoms();
  int len_subgraphs = num_subgraphs();

  for (int i = 0; i < len_atoms; i++) {
    if (atoms[i] == other && size() > 1) {
      paths.push_back({i + len_subgraphs});
    }
  }

  for (int i = 0; i < len_subgraphs; i++) {
    if (subgraphs[i].contains(other)) {
      auto r = subgraphs[i].get_paths_to(other);
      for (auto& v : r)
        v.insert(v.begin(), i);
      copy(r.begin(), r.end(), back_inserter(paths));
    }
  }

  return paths;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(
    const AEGraph& other) const {
  // returns all paths in the tree that lead to a subgraph like <other>
  std::vector<std::vector<int>> paths;
  int len_subgraphs = num_subgraphs();
  for (int i = 0; i < len_subgraphs; i++) {
    if (subgraphs[i] == other && size() > 1) {
      paths.push_back({i});
    } else {
      auto r = subgraphs[i].get_paths_to(other);
      for (auto& v : r)
        v.insert(v.begin(), i);
      copy(r.begin(), r.end(), back_inserter(paths));
    }
  }
  return paths;
}

std::vector<std::pair<std::string, int>> AEGraph::get_atoms(int level) {
  std::vector<std::pair<std::string, int>> level_atoms;
  int atoms_length = num_atoms();
  int subgraphs_length = num_subgraphs();
  for (int i = 0; i < atoms_length; i++) {
    level_atoms.push_back(std::make_pair(atoms[i], level));
  }
  for (int i = 0; i < subgraphs_length; i++) {
    std::vector<std::pair<std::string, int>> result;
    result = subgraphs[i].get_atoms(level + 1);
    level_atoms.insert(level_atoms.end(), result.begin(), result.end());
  }
  return level_atoms;
}
std::vector<std::vector<int>> AEGraph::possible_double_cuts() const {
  // functia este construita dupa modelul "get_paths_to", adica ne
  // deplasam prin graf si ne intoarcem recursiv pe caile care
  // ne intereseaza, adicam cele care contin "noduri" cu subgrafuri fara atomi,
  // dar cu un singur subgraf, deoarece din aceste noduri
  // se poate face double_cut
  std::vector<std::vector<int>> out_atoms;
  int subgraphs_length = num_subgraphs();

  for (int i = 0; i < subgraphs_length; i++) {
    if (subgraphs[i].num_atoms() == 0 && subgraphs[i].num_subgraphs() == 1) {
      out_atoms.push_back({i});
    }
  }
  for (int i = 0; i < subgraphs_length; i++) {
    auto r = subgraphs[i].possible_double_cuts();
    for (auto& v : r) {
      v.insert(v.begin(), i);
    }
    copy(r.begin(), r.end(), back_inserter(out_atoms));
  }

  return out_atoms;
}

AEGraph AEGraph::double_cut(std::vector<int> where) const {
  // tmp_graph_2 este un graf "temporar", cu ajutorul lui ne deplasam
  // prin subgrafurile de pe drumul dat ca parametru. Ne oprim cu un
  // pas inainte de final intrucat va trebui sa fim pozitionati
  // cu 2 nivele mai sus de locul de double_cut
  AEGraph tmp_graph = *this;
  AEGraph* tmp_graph_2 = &tmp_graph;
  for (unsigned int i = 0; i < where.size() - 1; i++) {
    tmp_graph_2 = &(tmp_graph_2->subgraphs[where[i]]);
  }
  int where_size = where.size() - 1;

  // reconstruim graful dupa procesul de double cut
  for (int i = 0;
       i <
       tmp_graph_2->subgraphs[where[where_size]].subgraphs[0].num_subgraphs();
       i++) {
    tmp_graph_2->subgraphs.push_back(
        tmp_graph_2->subgraphs[where[where_size]].subgraphs[0].subgraphs[i]);
  }
  for (int i = 0;
       i < tmp_graph_2->subgraphs[where[where_size]].subgraphs[0].num_atoms();
       i++) {
    tmp_graph_2->atoms.push_back(
        tmp_graph_2->subgraphs[where[where_size]].subgraphs[0].atoms[i]);
  }

  // stergem subgrafurile din locul unde se aflau inainte de double cut
  tmp_graph_2->subgraphs.erase(tmp_graph_2->subgraphs.begin() +
                               where[where_size]);
  return tmp_graph;
}

std::vector<std::vector<int>> AEGraph::possible_erasures(int level) const {
  // la fel, ca la possible double_cuts, in matricea r vom avea pe fiecare rand
  // un path catre un nod de pe nivel par.
  std::vector<std::vector<int>> path;
  int atoms_length = num_atoms();
  int subgraphs_length = num_subgraphs();
  if (atoms_length + subgraphs_length > 1 || level == -1) {
    for (int i = 0; i < subgraphs_length; i++) {
      if ((level + 1) % 2 == 0) {
        path.push_back({i});
      }
    }
    for (int i = 0; i < atoms_length; i++) {
      if ((level + 1) % 2 == 0) {
        path.push_back({i + subgraphs_length});
      }
    }
  }
  for (int i = 0; i < subgraphs_length; i++) {
    auto r = subgraphs[i].possible_erasures(level + 1);
    for (auto& v : r) {
      v.insert(v.begin(), i);
    }
    copy(r.begin(), r.end(), back_inserter(path));
  }

  return path;
}

AEGraph AEGraph::erase(std::vector<int> where) const {
  // cu graful temporar ne deplasam pana la parintele nodului
  // de la care trebuie sa facem erase, iar apoi stergem atomii
  // SAU subgrafurile acestui nod.
  AEGraph tmp_graph = *this;
  AEGraph* tmp_graph_2 = &tmp_graph;
  for (unsigned int i = 0; i < where.size() - 1; i++) {
    tmp_graph_2 = &(tmp_graph_2->subgraphs[where[i]]);
  }
  int where_size = where.size() - 1;
  if (where[where_size] < tmp_graph_2->num_subgraphs()) {
    tmp_graph_2->subgraphs.erase(tmp_graph_2->subgraphs.begin() +
                                 where[where_size]);
  } else {
    tmp_graph_2->atoms.erase(tmp_graph_2->atoms.begin() + where[where_size] -
                             tmp_graph_2->num_subgraphs());
  }
  return tmp_graph;
}
std::vector<std::vector<int>> AEGraph::get_all_paths() const {
  // returns all paths in the tree that lead to every atom
  std::vector<std::vector<int>> paths;

  int len_atoms = num_atoms();
  int len_subgraphs = num_subgraphs();

  for (int i = 0; i < len_atoms; i++) {
    paths.push_back({i + len_subgraphs});
  }

  for (int i = 0; i < len_subgraphs; i++) {
    auto r = subgraphs[i].get_all_paths();
    for (auto& v : r)
      v.insert(v.begin(), i);
    copy(r.begin(), r.end(), back_inserter(paths));
  }

  return paths;
}
bool AEGraph::is_equal(std::vector<int> vector_1,
                       std::vector<int> vector_2) const {
  if (vector_1.size() != vector_2.size()) {
    return false;
  }
  for (unsigned int i = 0; i < vector_1.size(); i++) {
    if (vector_1[i] != vector_2[i]) {
      return false;
    }
  }
  return true;
}

std::vector<std::vector<int>> AEGraph::possible_deiterations() const {
  // obtinem drumurile pana la atomi
  std::vector<std::vector<int>> path = this->get_all_paths();
  AEGraph tmp_graph = *this;
  AEGraph* tmp_graph_2 = &tmp_graph;
  std::vector<std::vector<int>> deiterate_path;
  // eliminam ultima muchie, intrucat dorim sa ne situam deaspura atomilor
  for (unsigned int i = 0; i < path.size(); i++) {
    path[i].pop_back();
    if (path[i].size() == 0) {
      path.erase(path.begin() + i);
      i--;
    }
  }
  // eliminam duplicatele
  for (unsigned int i = 0; i < path.size(); i++) {
    for (unsigned int j = i + 1; j < path.size(); j++) {
      if (is_equal(path[i], path[j])) {
        path.erase(path.begin() + j);
      }
    }
  }
  // pentru atomii si subgrafurile radacinei verificam existenta altor
  // drumuri spre atomi/subgrafuri identice
  int okay = -1;
  std::vector<std::vector<int>> atoms_path, subgraphs_path;
  for (int z = 0; z < tmp_graph_2->num_atoms(); z++) {
    atoms_path = get_paths_to(tmp_graph_2->atoms[z]);
    for (unsigned int ji = 0; ji < atoms_path.size(); ji++) {
      if (!is_equal({z + num_subgraphs()}, atoms_path[ji])) {
        if (1 < atoms_path[ji].size()) {
          deiterate_path.push_back(atoms_path[ji]);
          okay = z;
        }
      }
    }
  }
  for (int z = 0; z < tmp_graph_2->num_subgraphs(); z++) {
    subgraphs_path = get_paths_to(tmp_graph_2->subgraphs[z]);
    for (unsigned int ji = 0; ji < subgraphs_path.size(); ji++) {
      if (!is_equal({z}, subgraphs_path[ji])) {
        if (1 < subgraphs_path[ji].size()) {
          deiterate_path.push_back(subgraphs_path[ji]);
          okay = z;
        }
      }
    }
  }
  // realizam procesul de mai sus pentru fiecare nod din graf
  for (unsigned int i = 0; i < path.size(); i++) {
    tmp_graph_2 = &tmp_graph;
    int curr_path_size = path[i].size();
    std::vector<int> curr_path;

    for (int j = 0; j < curr_path_size; j++) {
      if (path[i][j] == okay) {
        continue;
      }
      tmp_graph_2 = &(tmp_graph_2->subgraphs[path[i][j]]);
      // folosim variabila okay pentru a selecta subgraful maximal
      int okay = 1;
      curr_path.push_back(path[i][j]);
      for (int z = 0; z < tmp_graph_2->num_atoms(); z++) {
        atoms_path = get_paths_to(tmp_graph_2->atoms[z]);
        curr_path.push_back(z + tmp_graph_2->num_subgraphs());
        for (unsigned int ji = 0; ji < atoms_path.size(); ji++) {
          if (!is_equal(curr_path, atoms_path[ji])) {
            if (curr_path.size() < atoms_path[ji].size()) {
              if (std::find(deiterate_path.begin(), deiterate_path.end(),
                            atoms_path[ji]) == deiterate_path.end()) {
                deiterate_path.push_back(atoms_path[ji]);
                okay = z;
              }
            }
          }
        }
        curr_path.pop_back();
      }
      for (int z = 0; z < tmp_graph_2->num_subgraphs(); z++) {
        subgraphs_path = get_paths_to(tmp_graph_2->subgraphs[z]);
        curr_path.push_back(z);
        for (unsigned int ji = 0; ji < subgraphs_path.size(); ji++) {
          if (!is_equal(curr_path, subgraphs_path[ji])) {
            if (curr_path.size() < subgraphs_path[ji].size()) {
              if (std::find(deiterate_path.begin(), deiterate_path.end(),
                            subgraphs_path[ji]) == deiterate_path.end()) {
                deiterate_path.push_back(subgraphs_path[ji]);
                okay = z;
              }
            }
          }
        }
        curr_path.pop_back();
      }
      if (okay >= 0) {
        return deiterate_path;
      }
    }
  }
  return deiterate_path;
}

AEGraph AEGraph::deiterate(std::vector<int> where) const {
  return erase(where);
}
#endif  // _AEGRAPH_H_
