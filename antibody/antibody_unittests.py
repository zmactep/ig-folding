#!/usr/bin/env python

import unittest
from antibody import *

class TestAntibody(unittest.TestCase):
    def test_numbering(self):
        with self.assertRaises(SystemExit) as cm:
            Numbering.process_region('AQAAAAAAAAABBBBLAC','VL','FR1')
            self.assertEqual(cm.exception.code, 1)
        with self.assertRaises(SystemExit) as cm:
            Numbering.process_region('AAAAAAAAAAAAAAAAAAAAQAAAAAAAAAABBBBLAC','VL','FR1')
            self.assertEqual(cm.exception.code, 1)
        with self.assertRaises(SystemExit) as cm:
            Numbering.process_region('AAAAAAAAAAAAAAAAAAAAQAAAAAAAAABBBBLAC','VL','FR1')
            self.assertEqual(cm.exception.code, 1)
        with self.assertRaises(SystemExit) as cm:
            Numbering.process_region('AAAAAAAAAAAAAAAAAAAAA','VL','FR1')
            self.assertEqual(cm.exception.code, 1)
        for s in ['AA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA']:
            for c in ['VL','VH']:
                for r in ['FR1','FR2','FR3','CDR1','CDR2','CDR3']:
                    with self.assertRaises(SystemExit) as cm:
                        Numbering.process_region(s,c,r)
                        self.assertEqual(cm.exception.code, 1)
        # VL FR1
        self.assertEqual(Numbering.process_region('AQAAAAAAAAAABBBBLAC','VL','FR1'),('AQAAAAAAAAAABBBBLAC','5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'))
        self.assertEqual(Numbering.process_region('AAQAAAAAAAAABBBBLAC','VL','FR1'),('AAQAAAAAAAAABBBBLAC','4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23'))

        self.assertEqual(Numbering.process_region('AAQAAAAAAAAAABBBBLAC','VL','FR1'),('AAQAAAAAAAAAABBBBLAC','4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'))
        self.assertEqual(Numbering.process_region('AAAQAAAAAAAAABBBBLAC','VL','FR1'),('AAAQAAAAAAAAABBBBLAC','3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23'))

        self.assertEqual(Numbering.process_region('AAAAAQAAAAAAAAAABBBBLAC','VL','FR1'),('AAAAAQAAAAAAAAAABBBBLAC'
            ,'1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'))
        self.assertEqual(Numbering.process_region('AAAAAQAAAAAAAAABBBBLAC','VL','FR1'),('AAAAAQAAAAAAAAABBBBLAC'
            ,'1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23'))
        self.assertEqual(Numbering.process_region('AAAAAAQAAAAAAAAAABBBBLAC','VL','FR1'),('AAAAAQAAAAAAAAAABBBBLAC'
            ,'1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'))
        self.assertEqual(Numbering.process_region('AAAAAAAQAAAAAAAAABBBBLAC','VL','FR1'),('AAAAAQAAAAAAAAABBBBLAC'
            ,'1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23'))
        # VL FR2
        self.assertEqual(Numbering.process_region('A'*15,'VL','FR2'),('A'*15,'35,36,37,38,39,40,41,42,43,44,45,46,47,48,49'))
        # VL FR3
        self.assertEqual(Numbering.process_region('A'*32,'VL','FR3'),
                ('A'*32,'57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88'))
        self.assertEqual(Numbering.process_region('A'*33,'VL','FR3'),
                ('A'*33,'57,58,59,60,61,62,63,64,65,66,66A,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88'))
        self.assertEqual(Numbering.process_region('A'*34,'VL','FR3'),
                ('A'*34,'57,58,59,60,61,62,63,64,65,66,66A,66B,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88'))
        # VL FR4
        self.assertEqual(Numbering.process_region('A'*12,'VL','FR4'),
                ('A'*12,'98,99,100,101,102,103,104,105,106,107,108,109'))
        # VL CDR1
        self.assertEqual(Numbering.process_region('A'*8,'VL','CDR1'),
                ('A'*8,'24,25,26,27,28,29,30,34'))
        self.assertEqual(Numbering.process_region('A'*9,'VL','CDR1'),
                ('A'*9,'24,25,26,27,28,29,30,33,34'))
        self.assertEqual(Numbering.process_region('A'*11,'VL','CDR1'),
                ('A'*11,'24,25,26,27,28,29,30,31,32,33,34'))
        self.assertEqual(Numbering.process_region('A'*12,'VL','CDR1'),
                ('A'*12,'24,25,26,27,28,29,30,30A,31,32,33,34'))
        self.assertEqual(Numbering.process_region('A'*15,'VL','CDR1'),
                ('A'*15,'24,25,26,27,28,29,30,30A,30B,30C,30D,31,32,33,34'))
        self.assertEqual(Numbering.process_region('A'*17,'VL','CDR1'),
                ('A'*17,'24,25,26,27,28,29,30,30A,30B,30C,30D,30E,30F,31,32,33,34'))
        # VL CDR2
        self.assertEqual(Numbering.process_region('A'*7,'VL','CDR2'),
                ('A'*7,'50,51,52,53,54,55,56'))
        self.assertEqual(Numbering.process_region('A'*8,'VL','CDR2'),
                ('A'*8,'50,51,52,53,54,54A,55,56'))
        self.assertEqual(Numbering.process_region('A'*9,'VL','CDR2'),
                ('A'*9,'50,51,52,53,54,54A,54B,55,56'))
        self.assertEqual(Numbering.process_region('A'*11,'VL','CDR2'),
                ('A'*11,'50,51,52,53,54,54A,54B,54C,54D,55,56'))
        # VL CDR3
        self.assertEqual(Numbering.process_region('A'*5,'VL','CDR3'),
                ('A'*5,'89,90,91,92,97'))
        self.assertEqual(Numbering.process_region('A'*6,'VL','CDR3'),
                ('A'*6,'89,90,91,92,93,97'))
        self.assertEqual(Numbering.process_region('A'*8,'VL','CDR3'),
                ('A'*8,'89,90,91,92,93,94,95,97'))
        self.assertEqual(Numbering.process_region('A'*9,'VL','CDR3'),
                ('A'*9,'89,90,91,92,93,94,95,96,97'))
        self.assertEqual(Numbering.process_region('A'*10,'VL','CDR3'),
                ('A'*10,'89,90,91,92,93,94,95,95A,96,97'))
        self.assertEqual(Numbering.process_region('A'*11,'VL','CDR3'),
                ('A'*11,'89,90,91,92,93,94,95,95A,95B,96,97'))
        self.assertEqual(Numbering.process_region('A'*14,'VL','CDR3'),
                ('A'*14,'89,90,91,92,93,94,95,95A,95B,95C,95D,95E,96,97'))
        self.assertEqual(Numbering.process_region('A'*15,'VL','CDR3'),
                ('A'*15,'89,90,91,92,93,94,95,95A,95B,95C,95D,95E,95F,96,97'))
        # VH FR1
        self.assertEqual(Numbering.process_region('A'*16,'VH','FR1'),
                ('A'*16,'10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25'))
        self.assertEqual(Numbering.process_region('A'*17,'VH','FR1'),
                ('A'*17,'9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25'))
        self.assertEqual(Numbering.process_region('A'*22,'VH','FR1'),
                ('A'*22,'4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25'))
        self.assertEqual(Numbering.process_region('A'*24,'VH','FR1'),
                ('A'*24,'2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25'))
        self.assertEqual(Numbering.process_region('A'*25,'VH','FR1'),
                ('A'*25,'1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25'))
        self.assertEqual(Numbering.process_region('A'*26,'VH','FR1'),
                ('A'*25,'1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25'))
        # VH FR2
        self.assertEqual(Numbering.process_region('A'*14,'VH','FR2'),
                ('A'*14,'36,37,38,39,40,41,42,43,44,45,46,47,48,49'))
        # VH FR3
        self.assertEqual(Numbering.process_region('A'*30,'VH','FR3'),
                ('A'*30,'66,67,68,69,70,71,72,73,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94'))
        self.assertEqual(Numbering.process_region('A'*31,'VH','FR3'),
                ('A'*31,'66,67,68,69,70,71,72,73,74,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94'))
        self.assertEqual(Numbering.process_region('A'*32,'VH','FR3'),
                ('A'*32,'66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94'))
        # VH FR4
        self.assertEqual(Numbering.process_region('A'*12,'VH','FR4'),
                ('A'*12,'103,104,105,106,107,108,109,110,111,112,113,114'))
        # VH CDR1
        self.assertEqual(Numbering.process_region('A'*6,'VH','CDR1'),
                ('A'*6,'26,27,32,33,34,35'))
        self.assertEqual(Numbering.process_region('A'*7,'VH','CDR1'),
                ('A'*7,'26,27,28,32,33,34,35'))
        self.assertEqual(Numbering.process_region('A'*9,'VH','CDR1'),
                ('A'*9,'26,27,28,29,30,32,33,34,35'))
        self.assertEqual(Numbering.process_region('A'*10,'VH','CDR1'),
                ('A'*10,'26,27,28,29,30,31,32,33,34,35'))
        self.assertEqual(Numbering.process_region('A'*11,'VH','CDR1'),
                ('A'*11,'26,27,28,29,30,31,31A,32,33,34,35'))
        self.assertEqual(Numbering.process_region('A'*13,'VH','CDR1'),
                ('A'*13,'26,27,28,29,30,31,31A,31B,31C,32,33,34,35'))
        # VH CDR2
        self.assertEqual(Numbering.process_region('A'*12,'VH','CDR2'),
                ('A'*12,'50,51,52,57,58,59,60,61,62,63,64,65'))
        self.assertEqual(Numbering.process_region('A'*13,'VH','CDR2'),
                ('A'*13,'50,51,52,56,57,58,59,60,61,62,63,64,65'))
        self.assertEqual(Numbering.process_region('A'*15,'VH','CDR2'),
                ('A'*15,'50,51,52,54,55,56,57,58,59,60,61,62,63,64,65'))
        self.assertEqual(Numbering.process_region('A'*16,'VH','CDR2'),
                ('A'*16,'50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65'))
        self.assertEqual(Numbering.process_region('A'*17,'VH','CDR2'),
                ('A'*17,'50,51,52,52A,53,54,55,56,57,58,59,60,61,62,63,64,65'))
        self.assertEqual(Numbering.process_region('A'*20,'VH','CDR2'),
                ('A'*20,'50,51,52,52A,52B,52C,52D,53,54,55,56,57,58,59,60,61,62,63,64,65'))
        self.assertEqual(Numbering.process_region('A'*22,'VH','CDR2'),
                ('A'*22,'50,51,52,52A,52B,52C,52D,52E,52F,53,54,55,56,57,58,59,60,61,62,63,64,65'))
        # VH CDR3
        self.assertEqual(Numbering.process_region('A'*3,'VH','CDR3'),
                ('A'*3,'95,96,97'))
        self.assertEqual(Numbering.process_region('A'*4,'VH','CDR3'),
                ('A'*4,'95,96,97,98'))
        self.assertEqual(Numbering.process_region('A'*6,'VH','CDR3'),
                ('A'*6,'95,96,97,98,99,100'))
        self.assertEqual(Numbering.process_region('A'*8,'VH','CDR3'),
                ('A'*8,'95,96,97,98,99,100,101,102'))
        self.assertEqual(Numbering.process_region('A'*9,'VH','CDR3'),
                ('A'*9,'95,96,97,98,99,100,100A,101,102'))
        self.assertEqual(Numbering.process_region('A'*18,'VH','CDR3'),
                ('A'*18,'95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,101,102'))
        self.assertEqual(Numbering.process_region('A'*34,'VH','CDR3'),
                ('A'*34,'95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,'+
                    '100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,100X,100Y,100Z,101,102'))


        def current_test(self):
            pass


if __name__ == '__main__':
    unittest.main()
