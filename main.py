import functions as func

# spacerpady controls the white space in-between the text fields
spacerpady = 20

# font_size controls the size font, text_font controlls the type of font display for the text fields
font_size = 10
text_font = "Segoe UI"

# NOTE Changing formatting options can lead to strange spacing as the .config() options for the individual text fields, labels, and buttons were designed with 
# these default values in mind, however I included them in this main code to be able to be changed to try to account for user preference

func.main(spacerpady, font_size, text_font)