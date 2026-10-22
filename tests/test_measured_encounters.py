import tempfile
import unittest
from pathlib import Path
import numpy as np
import openpyxl
from measured_encounters import batch_id, load_rows, weights, design, predict


class MeasuredEncounterTests(unittest.TestCase):
    def test_pooled_date_stays_one_group(self):
        self.assertEqual(batch_id('x/20220323_25_N'),'20220323_25')
        self.assertEqual(batch_id('x\\211020 image 4'),'20211020')
        with self.assertRaises(ValueError):batch_id('no experimental date')

    def test_molecule_and_batch_weights(self):
        rows=[{'batch':'a','molecule':'x'},{'batch':'a','molecule':'x'},
              {'batch':'a','molecule':'y'},{'batch':'b','molecule':'z'}]
        np.testing.assert_allclose(weights(rows),[.125,.125,.25,.5])

    def test_strict_boolean_and_missing_filters(self):
        w=openpyxl.Workbook();s=w.active;s.title='CTCF';s.append(['header']*15)
        good=[1,'x/211020','1.tif',True,'N',False,False,.05,100,5,7,30,40,'1-10',None]
        s.append(good);bad=good.copy();bad[3]=1;s.append(bad)
        bad=good.copy();bad[8]=None;s.append(bad)
        bad=good.copy();bad[4]='C';s.append(bad)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'test.xlsx';w.save(p);rows,counts=load_rows(p)
        self.assertEqual(len(rows),1);self.assertEqual(counts['CTCF']['not_boolean'],1)
        self.assertEqual(counts['CTCF']['missing_or_invalid_state'],1)
        self.assertNotIn('attempts',rows[0]);self.assertNotIn('residence',rows[0])

    def test_secondary_constant_uses_orientation(self):
        rows=[{'batch':str(i//2),'molecule':str(i),'side':i%2,'y':i%2} for i in range(8)]
        p=predict(rows,rows,'constant',secondary=True)
        self.assertTrue(np.all(p[::2]<.5));self.assertTrue(np.all(p[1::2]>.5))

    def test_features_and_probability_controls(self):
        rows=[{'batch':str(i//2),'molecule':str(i),'force':.01+.01*i,'elapsed':100*i,
               'loop':i+1,'side':0,'y':i%2} for i in range(8)]
        self.assertEqual(design(rows,'state_logistic',False).shape,(8,4))
        for model in ['constant','published_exponential','force_logistic','force_cubic','state_logistic']:
            p=predict(rows,rows,model)
            self.assertTrue(np.all(np.isfinite(p)));self.assertTrue(np.all((p>0)&(p<1)))

if __name__=='__main__':unittest.main()
