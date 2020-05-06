# -*- coding: utf-8 -*-
# Copyright (c) Ezcad Development Team. All Rights Reserved.

from qtpy.QtCore import Signal
from qtpy.QtWidgets import QLabel, QLineEdit, QGroupBox, QGridLayout, \
    QRadioButton, QButtonGroup, QHBoxLayout
from ezcad.config.base import _
from ezcad.widgets.ezdialog import EasyDialog


class Dialog(EasyDialog):
    NAME = _("Zoeppritz modeling")
    sig_start = Signal(tuple, str, str, str, str, str)

    def __init__(self, parent=None):
        EasyDialog.__init__(self, parent)
        self.setup_page()

    def setup_page(self):
        model_group = QGroupBox(_('Half-space elastic model'))
        lbl_model = QLabel(_("Model"))
        lbl_vp = QLabel(_("Vp"))
        lbl_vs = QLabel(_("Vs"))
        lbl_ro = QLabel(_("Density"))
        lbl_upp = QLabel(_("Upper space"))
        self.le_upp_vp = QLineEdit()
        self.le_upp_vs = QLineEdit()
        self.le_upp_ro = QLineEdit()
        lbl_low = QLabel(_("Lower space"))
        self.le_low_vp = QLineEdit()
        self.le_low_vs = QLineEdit()
        self.le_low_ro = QLineEdit()
        model_layout = QGridLayout()
        model_layout.addWidget(lbl_model, 0, 0)
        model_layout.addWidget(lbl_vp, 0, 1)
        model_layout.addWidget(lbl_vs, 0, 2)
        model_layout.addWidget(lbl_ro, 0, 3)
        model_layout.addWidget(lbl_upp, 1, 0)
        model_layout.addWidget(self.le_upp_vp, 1, 1)
        model_layout.addWidget(self.le_upp_vs, 1, 2)
        model_layout.addWidget(self.le_upp_ro, 1, 3)
        model_layout.addWidget(lbl_low, 2, 0)
        model_layout.addWidget(self.le_low_vp, 2, 1)
        model_layout.addWidget(self.le_low_vs, 2, 2)
        model_layout.addWidget(self.le_low_ro, 2, 3)
        model_group.setLayout(model_layout)
        self.layout.addWidget(model_group)

        text = _("Incident angles")
        self.angles = self.create_lineedit(text)
        self.layout.addWidget(self.angles)

        lbl_reflection = QLabel(_('Reflection'))
        self.rb_pp = QRadioButton(_('PP'))
        self.rb_ps = QRadioButton(_('PS'))
        bg_reflection = QButtonGroup()
        bg_reflection.addButton(self.rb_pp)
        bg_reflection.addButton(self.rb_ps)
        hbox = QHBoxLayout()
        hbox.addWidget(lbl_reflection)
        hbox.addWidget(self.rb_pp)
        hbox.addWidget(self.rb_ps)
        self.layout.addLayout(hbox)

        lbl_complex = QLabel(_('Complex'))
        self.rb_amp = QRadioButton(_('amplitude'))
        self.rb_pha = QRadioButton(_('phase'))
        bg_complex = QButtonGroup()
        bg_complex.addButton(self.rb_amp)
        bg_complex.addButton(self.rb_pha)
        hbox = QHBoxLayout()
        hbox.addWidget(lbl_complex)
        hbox.addWidget(self.rb_amp)
        hbox.addWidget(self.rb_pha)
        self.layout.addLayout(hbox)

        lbl_equation = QLabel(_('Equation'))
        self.rb_linear = QRadioButton(_('Linear'))
        self.rb_quadratic = QRadioButton(_('Quadratic'))
        self.rb_zoeppritz = QRadioButton(_('Zoeppritz'))
        bg_equation = QButtonGroup()
        bg_equation.addButton(self.rb_linear)
        bg_equation.addButton(self.rb_quadratic)
        bg_equation.addButton(self.rb_zoeppritz)
        hbox = QHBoxLayout()
        hbox.addWidget(lbl_equation)
        hbox.addWidget(self.rb_linear)
        hbox.addWidget(self.rb_quadratic)
        hbox.addWidget(self.rb_zoeppritz)
        self.layout.addLayout(hbox)

        text = _("New point")
        self.new_point = self.create_lineedit(text)
        self.layout.addWidget(self.new_point)

        action = self.create_action()
        self.layout.addWidget(action)

        self.rb_linear.toggled.connect(lambda: self.set_equation(self.rb_linear))
        self.rb_quadratic.toggled.connect(lambda: self.set_equation(self.rb_quadratic))
        self.rb_zoeppritz.toggled.connect(lambda: self.set_equation(self.rb_zoeppritz))
        self.rb_zoeppritz.setChecked(True)

        self.rb_pp.toggled.connect(lambda: self.set_reflection(self.rb_pp))
        self.rb_ps.toggled.connect(lambda: self.set_reflection(self.rb_ps))
        self.rb_pp.setChecked(True)

        self.rb_amp.toggled.connect(lambda: self.set_complex(self.rb_amp))
        self.rb_pha.toggled.connect(lambda: self.set_complex(self.rb_pha))
        self.rb_amp.setChecked(True)

    def set_equation(self, rb):
        if rb.isChecked():
            if rb.text() == _('Linear'):
                self.equation = 'linear'
            elif rb.text() == _('Quadratic'):
                self.equation = 'quadratic'
            elif rb.text() == _('Zoeppritz'):
                self.equation = 'zoeppritz'
            else:
                raise ValueError("Unknown value")

    def set_reflection(self, rb):
        if rb.isChecked():
            if rb.text() == _('PP'):
                self.reflection = 'PP'
            elif rb.text() == _('PS'):
                self.reflection = 'PS'
            else:
                raise ValueError("Unknown value")

    def set_complex(self, rb):
        if rb.isChecked():
            if rb.text() == _('amplitude'):
                self.complex = 'amplitude'
            elif rb.text() == _('phase'):
                self.complex = 'phase'
            else:
                raise ValueError("Unknown value")

    def apply(self):
        vp1 = float(self.le_upp_vp.text())
        vs1 = float(self.le_upp_vs.text())
        ro1 = float(self.le_upp_ro.text())
        vp2 = float(self.le_low_vp.text())
        vs2 = float(self.le_low_vs.text())
        ro2 = float(self.le_low_ro.text())
        model = (vp1, vs1, ro1, vp2, vs2, ro2)
        inc_angles = self.angles.edit.text()
        equation = self.equation
        reflection = self.reflection
        complexity = self.complex
        object_name = self.new_point.edit.text()
        self.sig_start.emit(model, inc_angles, equation, reflection,
                            complexity, object_name)


def main():
    from qtpy.QtWidgets import QApplication
    app = QApplication([])
    test = Dialog()
    test.sig_start.connect(print)
    test.show()
    app.exec_()


if __name__ == '__main__':
    main()
