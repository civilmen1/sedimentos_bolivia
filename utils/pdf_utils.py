from xhtml2pdf import pisa
import io

def create_pdf(html_content):
    """
    Converts HTML content to a PDF in memory.
    Returns the PDF data as bytes.
    """
    pdf_buffer = io.BytesIO()
    pisa_status = pisa.CreatePDF(html_content, dest=pdf_buffer)

    if pisa_status.err:
        return None

    return pdf_buffer.getvalue()
