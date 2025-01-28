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
        self.fit_ph1 = True
        self.plot_legend = False
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
        self.print_top_axis = True
        self.print_right_axis = True
        self.print_left_axis = True
        self.print_bottom_axis = True
        self.print_background = True
        self.print_light_mode = True
        self.print_dataset_colours = False
        self.print_standard_colours = True
        self.print_stacked_plot = False
        self.print_auto_scale = False
        self.print_stacked_plot_repeat_axes = False
        self.print_label = False
        self.print_spc_linewidth = 2
        self.print_axes_linewidth = 2
        self.print_ticks_font_size = 10
        self.print_label_font_size = 12
        self.print_nmr_spectrum_aspect_ratio = 'auto'
        self.print_hsqc_peak_aspect_ratio = 1.33
        self.print_hsqc_multiplet_aspect_ratio = 1.33
        self.local_baseline_correction = False

    def make_config(self):
        config = configparser.ConfigParser()
        auto_plot = 'yes' if self.auto_plot is True else 'no'
        keep_zoom = 'yes' if self.keep_zoom is True else 'no'
        plot_legend = 'yes' if self.plot_legend is True else 'no'
        print_top_axis = 'yes' if self.print_top_axis is True else 'no'
        print_right_axis = 'yes' if self.print_right_axis is True else 'no'
        print_left_axis = 'yes' if self.print_left_axis is True else 'no'
        print_bottom_axis = 'yes' if self.print_bottom_axis is True else 'no'
        print_background = 'yes' if self.print_background is True else 'no'
        print_light_mode = 'yes' if self.print_light_mode is True else 'no'
        print_dataset_colours = 'yes' if self.print_dataset_colours is True else 'no'
        print_standard_colours = 'yes' if self.print_standard_colours is True else 'no'
        print_stacked_plot = 'yes' if self.print_stacked_plot is True else 'no'
        print_stacked_plot_repeat_axes = 'yes' if self.print_stacked_plot_repeat_axes is True else 'no'
        print_auto_scale = 'yes' if self.print_auto_scale is True else 'no'
        fit_ph1 = 'no' if self.fit_ph1 is False else 'yes'
        print_label = 'yes' if self.print_label is True else 'no'
        local_baseline_correction = 'yes' if self.local_baseline_correction is True else 'no'
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
                          'neg_col22': self.neg_col22,
                          'plot_legend': plot_legend}
        config['System'] = {'current_directory': self.current_directory,
                            'local_baseline_correction': local_baseline_correction,
                            'fit_ph1': fit_ph1}
        config['Print'] = {'print_top_axis': print_top_axis,
                           'print_right_axis': print_right_axis,
                           'print_left_axis': print_left_axis,
                           'print_bottom_axis': print_bottom_axis,
                           'print_background': print_background,
                           'print_light_mode': print_light_mode,
                           'print_dataset_colours': print_dataset_colours,
                           'print_standard_colours': print_standard_colours,
                           'print_spc_linewidth': self.print_spc_linewidth,
                           'print_axes_linewidth': self.print_axes_linewidth,
                           'print_ticks_font_size': self.print_ticks_font_size,
                           'print_label_font_size': self.print_label_font_size,
                           'print_stacked_plot': print_stacked_plot,
                           'print_stacked_plot_repeat_axes': print_stacked_plot_repeat_axes,
                           'print_auto_scale': print_auto_scale,
                           'print_nmr_spectrum_aspect_ratio': str(self.print_nmr_spectrum_aspect_ratio),
                           'print_hsqc_peak_aspect_ratio': str(self.print_hsqc_peak_aspect_ratio),
                           'print_multiplet_peak_aspect_ratio': str(self.print_hsqc_multiplet_aspect_ratio),
                           'print_label': print_label
                           }
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
                local_var = local_var.replace("plotlegend", "plot_legend")
                local_var = local_var.replace("fontsize", "font_size")
                local_var = local_var.replace("phasereferencecolour", "phase_reference_colour")
                local_var = local_var.replace("poscol", "pos_col")
                local_var = local_var.replace("negcol", "neg_col")
                local_var = local_var.replace("localbaselinecorrection", "local_baseline_correction")
                self.set_values(local_var, config[k][l])

    def set_auto_plot(self, value):
        self.auto_plot = True if value == "yes" else False

    def set_current_directory(self, value):
        self.current_directory = value

    def set_keep_zoom(self, value):
        self.keep_zoom = True if value == "yes" else False

    def set_plot_legend(self, value):
        self.plot_legend = True if value == "yes" else False

    def set_font_size(self, value):
        self.font_size = float(value)

    def set_phase_reference_colour(self, value):
        self.phase_reference_colour = value

    def set_print_top_axis(self, value):
        self.print_top_axis = True if value == "yes" else False

    def set_print_right_axis(self, value):
        self.print_right_axis = True if value == "yes" else False

    def set_print_left_axis(self, value):
        self.print_left_axis = True if value == "yes" else False

    def set_print_bottom_axis(self, value):
        self.print_bottom_axis = True if value == "yes" else False

    def set_print_background(self, value):
        self.print_background = True if value == "yes" else False

    def set_print_light_mode(self, value):
        self.print_light_mode = True if value == "yes" else False

    def set_print_dataset_colours(self, value):
        self.print_dataset_colours = True if value == "yes" else False

    def set_print_standard_colours(self, value):
        self.print_standard_colours = True if value == "yes" else False

    def set_print_stacked_plot(self, value):
        self.print_stacked_plot = True if value == "yes" else False

    def set_print_stacked_plot_repeat_axes(self, value):
        self.print_stacked_plot_repeat_axes = True if value == "yes" else False

    def set_print_label(self, value):
        self.print_label = True if value == "yes" else False

    def set_print_auto_scale(self, value):
        self.print_auto_scale = True if value == "yes" else False

    def set_print_spc_linewidth(self, value):
        self.print_spc_linewidth = int(value)

    def set_print_axes_linewidth(self, value):
        self.print_axes_linewidth = int(value)

    def set_print_ticks_font_size(self, value):
        self.print_ticks_font_size = int(value)

    def set_print_label_font_size(self, value):
        self.print_label_font_size = int(value)

    def set_print_nmr_spectrum_aspect_ratio(self, value):
        if value == 'auto' or value == 'a4_landscape' or value == 'a4_portrait':
           self.print_nmr_spectrum_aspect_ratio = value
        else:
            self.print_nmr_spectrum_aspect_ratio = float(value)

    def set_print_hsqc_peak_aspect_ratio(self, value):
        if value == 'auto' or value == 'a4_landscape' or value == 'a4_portrait':
            self.print_hsqc_peak_aspect_ratio = value
        else:
            self.print_hsqc_peak_aspect_ratio = float(value)

    def set_print_multiplet_peak_aspect_ratio(self, value):
        if value == 'auto' or value == 'a4_landscape' or value == 'a4_portrait':
            self.print_multiplet_peak_aspect_ratio = value
        else:
            self.print_multiplet_peak_aspect_ratio = float(value)

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

    def set_fit_ph1(self, value):
        self.fit_ph1 = False if value == "no" else True

    def set_local_baseline_correction(self, value):
        self.local_baseline_correction = True if value == "yes" else False

    def set_values(self, key, value):
        m_name = "self.set_" + key
        eval(m_name)(value)
