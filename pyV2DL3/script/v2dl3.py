import logging


def main():
    logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
    logging.info(
        """
            Run V2DL3 for VEGAS (.stage5) or EventDisplay (.anasum):


            For VEGAS flags:

                v2dl3-vegas --help


            For Eventdisplay flags:

                pyV2DL3/script/v2dl3_for_Eventdisplay.py --help


            For more information, see V2DL3/README.md
            or
            https://github.com/VERITAS-Observatory/V2DL3
            """
    )


if __name__ == "__main__":
    main()
