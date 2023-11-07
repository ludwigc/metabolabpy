"""
NMR spectrum display parameters
"""

from metabolabpy.nmr import nmrConfig


class DispPars:

    def __init__(self):
        self.cf = nmrConfig.NmrConfig()
        self.cf.read_config()
        self.background_colour = []
        self.foreground_colour = []
        if self.cf.mode == 'light':
            self.pos_col = "RGB"  # (r,g,b) or colour string (e.g. 'b')
            self.pos_col_rgb = (self.cf.pos_col10, self.cf.pos_col11, self.cf.pos_col12)
            self.neg_col = "RGB"  # (r,g,b) or colour string (e.g. 'r')
            self.neg_col_rgb = (self.cf.neg_col10, self.cf.neg_col11, self.cf.neg_col12)
            self.background_colour = (255 / 255, 255 / 255, 255 / 255)
            self.foreground_colour = (0 / 255, 0 / 255, 0 / 255)
        else:
            self.pos_col = "RGB"  # (r,g,b) or colour string (e.g. 'b')
            self.pos_col_rgb = (self.cf.pos_col20, self.cf.pos_col21, self.cf.pos_col22)
            self.neg_col = "RGB"  # (r,g,b) or colour string (e.g. 'r')
            self.neg_col_rgb = (self.cf.neg_col20, self.cf.neg_col21, self.cf.neg_col22)
            self.background_colour = (42 / 255, 42 / 255, 42 / 255)
            self.foreground_colour = (255 / 255, 255 / 255, 255 / 255)

        self.ph_ref_col = self.cf.phase_reference_colour
        self.ph_ref_ds = 1
        self.ph_ref_exp = 1
        self.n_levels = 8  # integer1
        self.min_level = 0.01  # percentage of max
        self.max_level = 0.1  # percentage of max
        self.axis_type1 = 'ppm'  # 'ppm' / 'Hz'
        self.axis_type2 = 'ppm'  # 'ppm' / 'Hz'
        self.display_spc = False  # True / False (current spectrum is always plotted)
        self.spc_offset = 0.0  # percentage of max (1D only)
        self.spc_scale = 1.0  # fold change (1D only)
        self.x_label = '$^1$H'  # [axis_type] added during plot
        self.y_label = '$^{13}$C'  # [axis_type] added during plot
        #self.print_top_axis =
        self.spc_label = ""  # spectrum title if not empty
        self.int1 = 0.0
        self.int2 = 0.0
        self.int3 = 0.0
        self.print_colours = []
        self.print_neg_colours = []
        self.print_background_colour = []
        self.print_foreground_colour = []
        self.std_pos_col1 = []
        self.std_neg_col1 = []
        self.std_pos_col1 = []
        self.std_neg_col1 = []
        self.print_pos_col_rgb = self.pos_col_rgb
        self.print_neg_col_rgb = self.neg_col_rgb
        self.colours = {
            0: "RGB",
            1: "Red",
            2: "Green",
            3: "Blue",
            4: "Black",
            5: "Cyan",
            6: "Magenta",
            7: "Yellow",
            8: "Gray",
            9: "darkRed",
            10: "darkGreen",
            11: "darkBlue",
            12: "darkCyan",
            13: "darkMagenta",
            "RGB": 0,
            "Red": 1,
            "Green": 2,
            "Blue": 3,
            "Black": 4,
            "Cyan": 5,
            "Magenta": 6,
            "Yellow": 7,
            "Gray": 8,
            "darkRed": 9,
            "darkGreen": 10,
            "darkBlue": 11,
            "darkCyan": 12,
            "darkMagenta": 13
        }
        self.colours2 = {
            0: "Red",
            1: "Green",
            2: "Blue",
            3: "Black",
            4: "Cyan",
            5: "Magenta",
            6: "Yellow",
            7: "Gray",
            8: "darkRed",
            9: "darkGreen",
            10: "darkBlue",
            11: "darkCyan",
            12: "darkMagenta",
            "Red": 0,
            "Green": 1,
            "Blue": 2,
            "Black": 3,
            "Cyan": 4,
            "Magenta": 5,
            "Yellow": 6,
            "Gray": 7,
            "darkRed": 8,
            "darkGreen": 9,
            "darkBlue": 10,
            "darkCyan": 11,
            "darkMagenta": 12
        }
        self.axes = {
            0: "ppm",
            1: "Hz",
            "ppm": 0,
            "Hz": 1
        }
        self.false_true = {
            0: False,
            1: True,
            False: 0,
            True: 1
        }
        # end __init__

    def __str__(self):
        return "Display parameters"
        # end __str__

    def init_print_colours(self):
        self.cf.read_config()
        if self.cf.print_light_mode == False:
            self.int1 = 1.0
            self.int2 = 0.6
            self.int3 = 0.3
        else:
            self.int1 = 0.4
            self.int2 = 0.8
            self.int3 = 0.5

        int1 = self.int1
        int2 = self.int2
        int3 = self.int3
        light_neg_diff = 0.5
        dark_neg_diff = 0.2
        if self.cf.print_light_mode == False:
            self.print_colours = [(int1, int1, 0.0),
                                 (0.0, int1, int1),
                                 (int1, 0.0, int1),
                                 (int2, int2, int1),
                                 (int1, int2, int2),
                                 (int2, int1, int2),
                                 (int1, int1, int3),
                                 (int2, int3, int3),
                                 (int3, int2, int3),
                                 (int3, int2, int2),
                                 (int2, int2, int3),
                                 (int2, int3, int2)]

            int1 = self.int1 - dark_neg_diff
            int2 = self.int2 - dark_neg_diff
            int3 = self.int3 - dark_neg_diff
            int1 = max(int1, 0.0)
            int2 = max(int2, 0.0)
            int3 = max(int3, 0.0)
            self.print_neg_colours = [(int1, int1, 0.0),
                                 (0.0, int1, int1),
                                 (int1, 0.0, int1),
                                 (int2, int2, int1),
                                 (int1, int2, int2),
                                 (int2, int1, int2),
                                 (int1, int1, int3),
                                 (int2, int3, int3),
                                 (int3, int2, int3),
                                 (int3, int2, int2),
                                 (int2, int2, int3),
                                 (int2, int3, int2)]

        else:
            self.print_colours = [(0.0, 0.0, int1),
                                 (int1, 0.0, 0.0),
                                 (0.0, int1, 0.0),
                                 (0.0, int1, int1),
                                 (int1, int1, 0.0),
                                 (int1, 0.0, int1),
                                 (int3, int3, int2),
                                 (int2, int3, int3),
                                 (int3, int2, int3),
                                 (int3, int2, int2),
                                 (int2, int2, int3),
                                 (int2, int3, int2)]

            int1 = self.int1 + light_neg_diff
            int2 = self.int2 + light_neg_diff
            int3 = self.int3 + light_neg_diff
            int1 = min(int1, 1.0)
            int2 = min(int2, 1.0)
            int3 = min(int3, 1.0)
            self.print_neg_colours = [(light_neg_diff, light_neg_diff, int1),
                                 (int1, light_neg_diff, light_neg_diff),
                                 (light_neg_diff, int1, light_neg_diff),
                                 (light_neg_diff, int1, int1),
                                 (int1, int1, light_neg_diff),
                                 (int1, light_neg_diff, int1),
                                 (int3, int3, int2),
                                 (int2, int3, int3),
                                 (int3, int2, int3),
                                 (int3, int2, int2),
                                 (int2, int2, int3),
                                 (int2, int3, int2)]


        self.std_pos_col1 = (self.cf.pos_col10, self.cf.pos_col11, self.cf.pos_col12)
        self.std_neg_col1 = (self.cf.neg_col10, self.cf.neg_col11, self.cf.neg_col12)
        self.std_pos_col2 = (self.cf.pos_col20, self.cf.pos_col21, self.cf.pos_col22)
        self.std_neg_col2 = (self.cf.neg_col20, self.cf.neg_col21, self.cf.neg_col22)
        if self.cf.print_light_mode == False:
            self.print_pos_col_rgb = self.std_pos_col2
            self.print_neg_col_rgb = self.std_neg_col2
            self.print_background_colour = (42 / 255, 42 / 255, 42 / 255)
            self.print_foreground_colour = (255 / 255, 255 / 255, 255 / 255)
        else:
            self.print_pos_col_rgb = self.std_pos_col1
            self.print_neg_col_rgb = self.std_neg_col1
            self.print_background_colour = (255 / 255, 255 / 255, 255 / 255)
            self.print_foreground_colour = (0 / 255, 0 / 255, 0 / 255)

        # end init_print_colours

