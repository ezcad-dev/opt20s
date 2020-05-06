# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
How to run? pytest this.py or python -m pytest this.py
"""


from qtpy.QtCore import Qt
from gopoint.tool.zoep_dialog import Dialog


def test_zoep_dialog(qtbot):
    dialog = Dialog()
    dialog.sig_start.connect(print)
    dialog.show()
    assert dialog.isVisible()
    # assert dialog.windowTitle() == _("Zoeppritz modeling")

    ms_between_clicks = 1000  # msec

    qtbot.addWidget(dialog)
    qtbot.waitForWindowShown(dialog)
    qtbot.wait(ms_between_clicks)

    qtbot.keyClicks(dialog.le_upp_vp, '3.0')
    qtbot.wait(ms_between_clicks)
    qtbot.keyClicks(dialog.le_upp_vs, '1.5')
    qtbot.wait(ms_between_clicks)
    qtbot.keyClicks(dialog.le_upp_ro, '2.3')
    qtbot.wait(ms_between_clicks)
    qtbot.keyClicks(dialog.le_low_vp, '2.0')
    qtbot.wait(ms_between_clicks)
    qtbot.keyClicks(dialog.le_low_vs, '1.0')
    qtbot.wait(ms_between_clicks)
    qtbot.keyClicks(dialog.le_low_ro, '2.0')
    qtbot.wait(ms_between_clicks)
    qtbot.keyClicks(dialog.angles.edit, '0-60(1)')
    qtbot.wait(ms_between_clicks)

    # qtbot.mouseClick(dialog.rb_ps, Qt.LeftButton)  # not work
    dialog.rb_ps.setChecked(True)
    qtbot.wait(ms_between_clicks)

    # dialog.new_point.edit.setText('r_m1_ps_zoe_amp')
    qtbot.keyClicks(dialog.new_point.edit, 'r_m1_ps_zoe_amp')
    qtbot.wait(ms_between_clicks)

    qtbot.mouseClick(dialog.btn_apply, Qt.LeftButton)
    qtbot.wait(ms_between_clicks)

