import subprocess as sb
import os
import scipy.io
from pyformlang.finite_automaton import EpsilonNFA, State, Symbol


def run_intersection(a, b, w):
    TEST_DIR = os.path.dirname ( (os.path.abspath(__file__)) )
    INTERSECTION = os.path.join(os.path.dirname(TEST_DIR), "build/masked")
    NFA_DIV_DIR = os.path.join(os.path.dirname(TEST_DIR), "data/chrobelias/")
    args = f'2 {NFA_DIV_DIR}/{a}/1.txt {NFA_DIV_DIR}/{a}/2.txt {NFA_DIV_DIR}/{b}/1.txt {NFA_DIV_DIR}/{b}/2.txt \
    {NFA_DIV_DIR}/{a}/states.txt {NFA_DIV_DIR}/{b}/states.txt'
    sb.run(f'{INTERSECTION} {args} {TEST_DIR}/{w}.mtx', shell=True)
    return f'{TEST_DIR}/{w}.mtx'


def mtx_to_automaton(f):
    m = scipy.io.mmread(f, spmatrix=True)
    m_coo = m.tocoo()

    automaton = EpsilonNFA()
    automaton.add_start_state(State(0))
    automaton.add_final_state(State(0))

    for row, col, val in zip(m_coo.row, m_coo.col, m_coo.data):
        match val:
            case 1:
                automaton.add_transition(State(row), Symbol("1"), State(col))
            case 2:
                automaton.add_transition(State(row), Symbol("0"), State(col))
            case 3:
                automaton.add_transition(State(row), Symbol("0"), State(col))
                automaton.add_transition(State(row), Symbol("1"), State(col))

    return automaton


def load_automaton(name):
    TEST_DIR = os.path.dirname(os.path.abspath(__file__))
    NFA_DIV_DIR = os.path.join(os.path.dirname(TEST_DIR), "data/chrobelias/")
    d = f"{NFA_DIV_DIR}/{name}/"

    m0 = scipy.io.mmread(f"{d}/1.txt", spmatrix=True).tocoo()
    m1 = scipy.io.mmread(f"{d}/2.txt", spmatrix=True).tocoo()
    st = scipy.io.mmread(f"{d}/states.txt", spmatrix=True).tocoo()

    automaton = EpsilonNFA()
    for r, c in zip(m0.row, m0.col):
        automaton.add_transition(State(r), Symbol("0"), State(c))
    for r, c in zip(m1.row, m1.col):
        automaton.add_transition(State(r), Symbol("1"), State(c))
    for s in set(st.row):
        automaton.add_start_state(State(s))
        automaton.add_final_state(State(s))

    return automaton


class TestIntersection_div3_div9:
    def setup_method(self):
        automaton_path = run_intersection("nfa-div-3", "nfa-div-9", "new-nfa-div-9")
        self.automaton = mtx_to_automaton(automaton_path)
        self.reference_automaton = load_automaton("nfa-div-3").get_intersection(load_automaton("nfa-div-9"))

    def test_equivalence(self):
        assert self.automaton.is_equivalent_to(self.reference_automaton)


class TestIntersection_div9_div15:
    def setup_method(self):
        automaton_path = run_intersection("nfa-div-9", "nfa-div-15", "new-nfa-div-45")
        self.automaton = mtx_to_automaton(automaton_path)
        self.reference_automaton = load_automaton("nfa-div-9").get_intersection(load_automaton("nfa-div-15"))

    def test_equivalence(self):
        assert self.automaton.is_equivalent_to(self.reference_automaton)


class TestIntersection_div3_div7:
    def setup_method(self):
        automaton_path = run_intersection("nfa-div-3", "nfa-div-7", "new-nfa-div-21")
        self.automaton = mtx_to_automaton(automaton_path)
        self.reference_automaton = load_automaton("nfa-div-3").get_intersection(load_automaton("nfa-div-7"))

    def test_equivalence(self):
        assert self.automaton.is_equivalent_to(self.reference_automaton)
