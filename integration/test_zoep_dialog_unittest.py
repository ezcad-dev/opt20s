# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.

import unittest
from qtpy.QtWidgets import QApplication
from gopoint.tool.zoep_dialog import Dialog


class Test(unittest.TestCase):
    def test_1(self):
        # Two half spaces elastic model
        vp1, vp2 = 3.0, 2.0
        vs1, vs2 = 1.5, 1.0
        ro1, ro2 = 2.3, 2.0

        app = QApplication([])
        dialog = Dialog()
        dialog.sig_start.connect(print)
        dialog.show()
        # app.exec_()  # blocks thread for GUI user input

        dialog.le_upp_vp.setText(str(vp1))
        dialog.le_upp_vs.setText(str(vs1))
        dialog.le_upp_ro.setText(str(ro1))
        dialog.le_low_vp.setText(str(vp2))
        dialog.le_low_vs.setText(str(vs2))
        dialog.le_low_ro.setText(str(ro2))
        dialog.angles.edit.setText('1,2,3')
        dialog.rb_ps.setChecked(True)
        dialog.rb_zoeppritz.setChecked(True)
        dialog.rb_amp.setChecked(True)
        dialog.new_point.edit.setText('r_m1_ps_zoe_amp')

        dialog.apply()
        app.exec_()  # blocks thread for GUI user input
        # The limitation is cannot simulate user typing after dialog pop-up.


if __name__ == '__main__':
    unittest.main()
