import os
import configparser


class NmrConfig:

    def __init__(self):
        self.auto_plot = False
        self.keep_zoom = True
        self.font_size = 13
        self.phase_reference_colour = 'Red'
        self.f_name = '.pyMetaboLab.config'
        self.home_dir = os.path.expanduser('~')
        self.config_file = os.path.join(self.home_dir, self.f_name)
        self.value = ""
        self.pos_col10 = 0.0
        self.pos_col11 = 0.0
        self.pos_col12 = 1.0
        self.neg_col10 = 1.0
        self.neg_col11 = 0.0
        self.neg_col12 = 0.0
        self.pos_col20 = 0.8
        self.pos_col21 = 0.8
        self.pos_col22 = 1.0
        self.neg_col20 = 1.0
        self.neg_col21 = 0.8
        self.neg_col22 = 0.8
        self.mode = 'system'
        self.current_directory = ''

    def make_config(self):
        config = configparser.ConfigParser()
        auto_plot = 'yes' if self.auto_plot is True else 'no'
        keep_zoom = 'yes' if self.keep_zoom is True else 'no'
        config['GUI'] = {'auto_plot': auto_plot,
                         'keep_zoom': keep_zoom,
                         'font_size': str(self.font_size),
                         'mode': self.mode}
        config['Disp'] = {'phase_reference_colour': self.phase_reference_colour,
                          'pos_col10': self.pos_col10,
                          'pos_col11': self.pos_col11,
                          'pos_col12': self.pos_col12,
                          'neg_col10': self.neg_col10,
                          'neg_col11': self.neg_col11,
                          'neg_col12': self.neg_col12,
                          'pos_col20': self.pos_col20,
                          'pos_col21': self.pos_col21,
                          'pos_col22': self.pos_col22,
                          'neg_col20': self.neg_col20,
                          'neg_col21': self.neg_col21,
                          'neg_col22': self.neg_col22}
        config['System'] = {'current_directory': self.current_directory}
        return config

    def save_config(self):
        config = self.make_config()
        with open(self.config_file, 'w') as config_file:
            config.write(config_file)

    def read_config(self):
        self.home_dir = os.path.expanduser('~')
        self.config_file = os.path.join(self.home_dir, self.f_name)
        config = configparser.ConfigParser()
        config.read(self.config_file)
        if len(config.sections()) == 0:
            config = self.make_config()
            self.save_config()

        for k in config.sections():
            for l in config[k]:
                local_var = l
                local_var = local_var.replace("autoplot", "auto_plot")
                local_var = local_var.replace("keepzoom", "keep_zoom")
                local_var = local_var.replace("fontsize", "font_size")
                local_var = local_var.replace("phasereferencecolour", "phase_reference_colour")
                local_var = local_var.replace("poscol", "pos_col")
                local_var = local_var.replace("negcol", "neg_col")
                self.set_values(local_var, config[k][l])

    def set_auto_plot(self, value):
        self.auto_plot = True if value == "yes" else False

    def set_current_directory(self, value):
        self.current_directory = value

    def set_keep_zoom(self, value):
        self.keep_zoom = True if value == "yes" else False

    def set_font_size(self, value):
        self.font_size = float(value)

    def set_phase_reference_colour(self, value):
        self.phase_reference_colour = value

    def set_pos_col10(self, value):
        self.pos_col10 = float(value)

    def set_pos_col11(self, value):
        self.pos_col11 = float(value)

    def set_pos_col12(self, value):
        self.pos_col12 = float(value)

    def set_neg_col10(self, value):
        self.neg_col10 = float(value)

    def set_neg_col11(self, value):
        self.neg_col11 = float(value)

    def set_neg_col12(self, value):
        self.neg_col12 = float(value)

    def set_pos_col20(self, value):
        self.pos_col20 = float(value)

    def set_pos_col21(self, value):
        self.pos_col21 = float(value)

    def set_pos_col22(self, value):
        self.pos_col22 = float(value)

    def set_neg_col20(self, value):
        self.neg_col20 = float(value)

    def set_neg_col21(self, value):
        self.neg_col21 = float(value)

    def set_neg_col22(self, value):
        self.neg_col22 = float(value)

    def set_mode(self, value):
        self.mode = value

    def set_values(self, key, value):
        m_name = "self.set_" + key
        eval(m_name)(value)
