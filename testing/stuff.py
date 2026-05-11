import os, sys

ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')
)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


from class1 import c1_matching_comparison as match

<<<<<<< HEAD
match.plot
=======
match.plot_sensitivity_study([r'outputs\Matching_concepts\Taildragger Boosted Piston_A_results.csv'], 
                             [r'outputs\Matching_concepts\Taildragger Boosted Piston_CL_results.csv'],
                             r'outputs')
>>>>>>> 19acaa9 (performing sensitivity study)
