
from wtforms import Form, StringField
from flask import Flask, request, render_template

from feynman import make_amplitude, ParseException, calculate
import traceback

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index(result="", report_visible=False, error=""):
    if request.method == "POST":
        try:
            result = calculate(request.form.get("config_str"), request.form.get("internal_momenta"))
        #except ParseException as e:
        #    error = e.text
        except Exception as e:
            traceback.print_exc()
        report_visible = True
    return render_template("calculate.html", result=result, report_visible=report_visible, error=error)
