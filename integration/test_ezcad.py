# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.
"""
GUI testing with pytest-qt
How to run? pytest this.py or python -m pytest this.py
"""

# import os
import copy
from datetime import datetime
from argparse import Namespace
from dotenv import load_dotenv
from qtpy.QtCore import Qt
from ezcad.app.mainwindow import MainWindow


def test_ezcad(qtbot):
    load_dotenv()
    # print('Load javaseis module', os.getenv('data_io/javaseis'))

    path = 'C:\\Users\\xinfa\\Documents\\code\\ezcad-dev\\'
    options = Namespace(
        plugins_paths=path+'plugins\\ezcad_plugins',
        working_directory=None,
        debug=False
    )

    ms_between_clicks = 1000  # msec

    window = MainWindow(options)
    window.setup()
    window.show()
    window.post_visible_setup()
    fn = "C:\\Users\\xinfa\\Documents\\test\\test{}.ezd".format(
        datetime.now().strftime("%Y%m%d%H%M%S"))
    window.create_new_project(fn)
    # qtbot.addWidget(window)
    # qtbot.wait(1000)

    idx = window.menubar_plus.toolbar.cb_mode.findText("Point Menubar")
    window.menubar_plus.toolbar.cb_mode.setCurrentIndex(idx)
    qtbot.wait(ms_between_clicks)

    pm = window.plugins_manager.getPluginByName("Point Menubar",
                                                category="PluginMenuBar")
    # print('test1', pm.plugin_object)
    # print('test2', pm.plugin_object.bartender)
    # qtbot click tool > then drop-down menu?
    dialog = pm.plugin_object.bartender.zoep_modeling()

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

    name1, name2 = 'r_m1_ps_zoe_amp', 'r_m2_ps_zoe_amp'
    # dialog.new_point.edit.setText('r_m1_ps_zoe_amp')
    qtbot.keyClicks(dialog.new_point.edit, name1)
    qtbot.wait(ms_between_clicks)

    qtbot.mouseClick(dialog.btn_apply, Qt.LeftButton)
    qtbot.wait(2*ms_between_clicks)

    # Change density and create new point
    dialog.le_low_ro.clear()
    qtbot.wait(ms_between_clicks)
    qtbot.keyClicks(dialog.le_low_ro, '2.1')
    qtbot.wait(ms_between_clicks)
    dialog.new_point.edit.clear()
    qtbot.wait(ms_between_clicks)
    qtbot.keyClicks(dialog.new_point.edit, name2)
    qtbot.wait(ms_between_clicks)
    qtbot.mouseClick(dialog.btn_ok, Qt.LeftButton)
    qtbot.wait(ms_between_clicks)

    # open plot viewer
    window.viewer.tabs.new_vs_plot()
    qtbot.wait(ms_between_clicks)
    window.current_viewer.apply_aspect('Auto')
    qtbot.wait(ms_between_clicks)

    # check object in tree - turn on in viewer
    window.treebase.object_items[name1].setCheckState(0, Qt.Checked)
    qtbot.wait(2*ms_between_clicks)
    window.treebase.object_items[name2].setCheckState(0, Qt.Checked)
    qtbot.wait(ms_between_clicks)

    # style editor, change color
    dob = window.database[name2]
    style = copy.deepcopy(dob.atom_style)
    style['face_color'] = (0, 0, 255, 255)
    dob.set_atom_style(atom_style=style, update_plots=True)

    # flash it - turn on/off
    window.treebase.object_items[name2].setCheckState(0, Qt.Unchecked)
    qtbot.wait(2*ms_between_clicks)
    window.treebase.object_items[name2].setCheckState(0, Qt.Checked)
    qtbot.wait(2*ms_between_clicks)
    window.treebase.object_items[name2].setCheckState(0, Qt.Unchecked)
    qtbot.wait(2*ms_between_clicks)
    window.treebase.object_items[name2].setCheckState(0, Qt.Checked)
    qtbot.wait(2*ms_between_clicks)
