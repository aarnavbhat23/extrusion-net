import itertools
import unittest
import numpy as np
from trajectory_reversibility import track_scores,reversal_test


class ReversibilityTests(unittest.TestCase):
    def test_reversal_and_gaps(self):
        t=np.array([0,1,3,4,8]);d=np.array([.1,.8,.4,.9,np.nan])
        score,count=track_scores(t,d,[.3,.6]);reverse,_=track_scores(8-t,d,[.3,.6])
        np.testing.assert_allclose(score,-reverse);self.assertEqual(count,[2,1,1])

    def test_duplicate_frames_rejected(self):
        with self.assertRaises(ValueError):track_scores([0,0],[.1,.2],[.3,.6])

    def test_reversible_orbit_calibration(self):
        x=np.array([[.3,.2],[.1,-.2],[-.4,.1],[.2,.3]])
        p=[reversal_test(x*np.array(s)[:,None])['adjusted_p'] for s in itertools.product((-1,1),repeat=4)]
        for alpha in (.05,.125,.25,.5):self.assertLessEqual(np.mean(np.array(p)<=alpha),alpha)

    def test_visible_and_hidden_cycle(self):
        t=np.arange(600);d=np.tile([.1,.5,.9],200)
        score,_=track_scores(t,d,[.3,.6]);self.assertGreater(score[0],.99)
        result=reversal_test(np.tile(score,(8,1)));self.assertTrue(result['gate_passed'])
        invisible,_=track_scores(t,np.ones(600),[.3,.6]);np.testing.assert_array_equal(invisible,0)
        self.assertFalse(reversal_test(np.tile(invisible,(8,1)))['gate_passed'])

    def test_insufficient_clusters(self):
        result=reversal_test(np.ones((4,3)));self.assertFalse(result['sufficient_clusters'])

    def test_stationary_reversible_chain_exact_flux(self):
        # Symmetric transition kernel: stationary uniform distribution satisfies detailed balance.
        p=np.array([[.7,.2,.1],[.2,.5,.3],[.1,.3,.6]])
        cycle=np.array([[0,1,-1],[-1,0,1],[1,-1,0]])
        for lag in (1,2,4):
            joint=np.linalg.matrix_power(p,lag)/3
            self.assertAlmostEqual(float(np.sum(joint*cycle)),0.,places=14)

    def test_duplicate_tracks_do_not_increase_cluster_weight(self):
        import csv
        import tempfile
        from pathlib import Path
        from study_trajectory_reversibility import analyze
        rows=[]
        for day in range(1,22):
            for track in (0,1):
                for frame in range(21):
                    rows.append({'movie_index':f'2020_01_{day:02d}_movie1','id':str(track),'t':frame,
                                 'dist':(.1,.5,.9)[(frame+track)%3]})
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'tracks.tsv'
            def write(data):
                with path.open('w') as handle:
                    writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(data)
            write(rows);original=analyze(path)
            write(rows+[{**r,'id':str(int(r['id'])+2)} for r in rows]);duplicated=analyze(path)
        self.assertEqual(original['test'],duplicated['test'])
        self.assertEqual(original['groups'],duplicated['groups'])
