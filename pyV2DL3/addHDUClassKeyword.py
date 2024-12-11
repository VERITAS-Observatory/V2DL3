from pyV2DL3.constant import HDUCLASS, HDUDOC, HDUVERS


def addHDUClassKeyword(hdu, class1, class2=None, class3=None, class4=None):
    hdu.header.set("HDUCLASS", HDUCLASS, "This FITS file follows the GADF data format")
    hdu.header.set("HDUDOC", HDUDOC)
    hdu.header.set("HDUVERS", HDUVERS, "DL3 specification version")
    hdu.header.set("HDUCLAS1", class1, "Primary extension class")
    if class2 is not None:
        hdu.header.set("HDUCLAS2", class2)
    if class3 is not None:
        hdu.header.set("HDUCLAS3", class3)
    if class4 is not None:
        hdu.header.set("HDUCLAS4", class4)
    return hdu
