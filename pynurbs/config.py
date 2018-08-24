class Settings(object):
    """
    Settings.

    :var float gtol: Geometric tolerance (default=1.0e-7).
    :var float ptol: Parametric tolerance (default=1.0e-9).
    :var float atol: Angular tolerance (default=1.0e-12).
    :var float ftol: Used for curve/surface flatness criteria in
        subdivision methods (default=1.0e-3).
    :var float mtol: Mesh tolerance (default=0.005).
    :var float part_tol: Part tolerance (default=0.005).
    :var float part_angle: Part angle limit in degrees (default=30.0).
    :var bool use_fortran: Option to use Fortran methods if available
        (default= *True*).
    var bool warnings: Option to print warning messages (default=*False*).
    """
    # Class variables for settings.
    gtol = 1.0e-7
    ptol = 1.0e-9
    atol = 1.0e-12
    ftol = 1.0e-3
    mtol = 0.005
    part_tol = 0.005
    part_angle = 30.0
    use_fortran = True
    warnings = False

    # @classmethod
    # def set_fortran_config(cls):
    #     """
    #     Set Fortran configuration values.
    #     """
    #     lib_config.gtol = cls.gtol
    #     lib_config.ptol = cls.ptol
    #     lib_config.atol = cls.atol
    #     lib_config.stol = cls.stol
    #     lib_config.ftol = cls.ftol
    #     lib_config.warnings = cls.warnings

    @classmethod
    def set_gtol(cls, gtol=1.0e-7):
        """
        Set the default geometric tolerance.

        :param float gtol: Geometric tolerance.
        """
        cls.gtol = float(gtol)
        # lib_config.gtol = cls.gtol

    @classmethod
    def set_ptol(cls, ptol=1.0e-12):
        """
        Set the default parametric tolerance.

        :param float ptol: Parametric tolerance.
        """
        cls.ptol = float(ptol)
        # lib_config.ptol = cls.ptol

    @classmethod
    def set_atol(cls, atol=1.0e-12):
        """
        Set the default angular tolerance.

        :param float atol: Angular tolerance.
        """
        cls.atol = float(atol)
        # lib_config.atol = cls.atol

    @classmethod
    def set_ftol(cls, ftol=1.0e-3):
        """
        Set the default tolerance for curve/surface flatness criteria.

        :param float ftol: Flatness tolerance.
        """
        cls.ftol = float(ftol)
        # lib_config.ftol = cls.ftol

    @classmethod
    def set_mtol(cls, mtol=0.005):
        """
        Set the default mesh tolerance.

        :param float mtol: Mesh tolerance.
        """
        cls.mtol = float(mtol)

    @classmethod
    def set_part_tol(cls, part_tol=0.005):
        """
        Set the default part tolerance.

        :param float part_tol: Part tolerance.
        """
        cls.part_tol = float(part_tol)

    @classmethod
    def set_part_angle(cls, part_angle=30.0):
        """
        Set the default part angle limit.

        :param float part_angle: Part angle limit.
        """
        cls.part_angle = float(part_angle)

    @classmethod
    def activate_fortran(cls):
        """
        Activate the use of Fortran methods.
        """
        cls.use_fortran = True

    @classmethod
    def deactivate_fortran(cls):
        """
        Deactivate the use of Fortran methods.
        """
        cls.use_fortran = False

    @classmethod
    def activate_warnings(cls):
        """
        Activate printing warning messages.
        """
        cls.warnings = True
        # lib_config.warnings = True

    @classmethod
    def deactivate_warnings(cls):
        """
        Deactivate printing warning messages.
        """
        cls.warnings = False
        # lib_config.warnings = False

# Set Fortran config variables.
# Settings.set_fortran_config()
