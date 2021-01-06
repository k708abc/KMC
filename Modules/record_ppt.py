from pptx import Presentation
from pptx.util import Inches, Pt
import os
import datetime


def rec_ppt(params, minute, second, img_names, hist_names, time, coverage, dir_name):
    ppt_name = dir_name + "Layer_analysis_results.pptx"
    if os.path.exists(ppt_name):
        prs = Presentation(ppt_name)
    else:
        prs = Presentation()
        prs.slide_height = Inches(7.5)
        prs.slide_width = Inches(13.33)
        title_slide_layout = prs.slide_layouts[0]
        slide = prs.slides.add_slide(title_slide_layout)
        title = slide.shapes.title
        title.text = "Layer analysis calculation results"
    # First page
    blank_slide_layout = prs.slide_layouts[6]
    slide = prs.slides.add_slide(blank_slide_layout)
    width = height = Inches(1)
    top = Inches(0)
    left = Inches(0.5)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    p = tf.add_paragraph()
    p.text = "Parameters"
    p.font.size = Pt(28)
    shapes = slide.shapes
    left = Inches(0.1)
    top = Inches(1.5)
    rows = 2
    cols = 9
    width = Inches(13)
    height = Inches(1)
    # record parameters
    table0 = shapes.add_table(rows, cols, left, top, width, height).table
    table0.cell(0, 0).text = "Num. cells"
    table0.cell(0, 1).text = "Z units"
    table0.cell(0, 2).text = "T (K)"
    table0.cell(0, 3).text = "kbT"
    table0.cell(0, 4).text = "dep. rate (ML/min)"
    table0.cell(0, 5).text = "dep. time (min)"
    table0.cell(0, 6).text = "prefactor"
    table0.cell(0, 7).text = "Cal. time (s)"
    table0.cell(0, 8).text = "Finishing time"
    table0.cell(1, 0).text = str(params.n_cell_init)
    table0.cell(1, 1).text = str(params.z_unit_init)
    table0.cell(1, 2).text = str(params.temperature)
    table0.cell(1, 3).text = str("{:.3g}".format(params.temperature_eV))
    table0.cell(1, 4).text = str(params.dep_rate)
    table0.cell(1, 5).text = str(params.dep_time)
    table0.cell(1, 6).text = str("{:.1E}".format(params.prefactor))
    table0.cell(1, 7).text = str(minute) + " m " + str(second) + " s"
    dt_now = datetime.datetime.now()
    table0.cell(1, 8).text = str(dt_now.strftime("%Y/%m/%d/ %H:%M:%S"))
    #
    left = Inches(0.1)
    top = Inches(3.5)
    rows = 3
    cols = 6
    width = Inches(13)
    height = Inches(1)

    table = shapes.add_table(rows, cols, left, top, width, height).table
    table.cell(0, 0).text = "Transformation"
    table.cell(0, 1).text = "Keep defects"
    table.cell(0, 2).text = "Put at first"
    table.cell(0, 3).text = "Cut event"
    table.cell(0, 4).text = "Rate limit"
    table.cell(0, 5).text = "Method"
    table.cell(1, 0).text = str(params.trans_check)
    table.cell(1, 1).text = str(params.keep_defect_check)
    table.cell(1, 2).text = str(params.first_put_check)
    table.cell(1, 3).text = str(params.cut_check)
    table.cell(1, 4).text = str(params.limit_check)
    table.cell(1, 5).text = str(params.method)
    table.cell(2, 0).text = str(params.transformation)
    table.cell(2, 1).text = str(params.num_defect)
    table.cell(2, 2).text = str(params.put_first)
    table.cell(2, 3).text = str(params.cut_number)
    table.cell(2, 4).text = str(params.limit_val)
    table.cell(2, 5).text = str("")
    #
    left = Inches(0.1)
    top = Inches(5)
    rows = 2
    cols = 12
    width = Inches(13)
    height = Inches(1)
    table = shapes.add_table(rows, cols, left, top, width, height).table
    table.cell(0, 1).text = "Ag base"
    table.cell(0, 2).text = "Ag-Si"
    table.cell(0, 3).text = "Si base"
    table.cell(0, 4).text = "Si(0-1)"
    table.cell(0, 5).text = "Si(1-2)"
    table.cell(0, 6).text = "Si(2-3)"
    table.cell(0, 7).text = "Si(3-4)"
    table.cell(0, 8).text = "Si(4-5)"
    table.cell(0, 9).text = "Si(inter)"
    table.cell(0, 10).text = "Si(intra)"
    table.cell(0, 11).text = "ES"
    # table.cell(0, 9).text = "Ag(top)"
    # table.cell(0, 10).text = "Trans."
    table.cell(1, 0).text = "Energy"
    table.cell(1, 1).text = str(params.binding_energies["Ag base"])
    table.cell(1, 2).text = str(params.binding_energies["AgSi"])
    table.cell(1, 3).text = str(params.binding_energies["Si base"])
    table.cell(1, 4).text = str(params.binding_energies["Si01"])
    table.cell(1, 5).text = str(params.binding_energies["Si12"])
    table.cell(1, 6).text = str(params.binding_energies["Si23"])
    table.cell(1, 7).text = str(params.binding_energies["Si34"])
    table.cell(1, 8).text = str(params.binding_energies["Si45"])
    table.cell(1, 9).text = str(params.binding_energies["Si_inter"])
    table.cell(1, 10).text = str(params.binding_energies["Si_intra"])
    table.cell(1, 11).text = str(params.binding_energies["ES"])
    # table.cell(1, 9).text = str(params.binding_energies["Agtop"])
    # table.cell(1, 10).text = str(params.transformation)
    #
    """
    left = Inches(0.5)
    top = Inches(6.5)
    rows = 2
    cols = 3
    width = Inches(8)
    height = Inches(0.5)

    table = shapes.add_table(rows, cols, left, top, width, height).table
    table.cell(0, 0).text = "Total energy"
    table.cell(0, 1).text = "1 ML"
    table.cell(0, 2).text = "Final"

    table.cell(1, 0).text = str('{:.3g}'.format(E))
    table.cell(1, 1).text = str('{:.3g}'.format(ML_check))
    table.cell(1, 2).text = str('{:.3g}'.format(final_check))
    """
    #
    width = height = Inches(1)
    top = Inches(6)
    left = Inches(0.5)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    p = tf.add_paragraph()
    p.text = "Comment: " + "\n" + params.comments
    p.font.size = Pt(20)

    #
    num_ims = len(img_names)
    imn = 0
    while imn < num_ims:
        slide = prs.slides.add_slide(blank_slide_layout)
        width = height = Inches(1)
        top = Inches(-0.1)
        left = Inches(0.5)
        txBox = slide.shapes.add_textbox(left, top, width, height)
        tf = txBox.text_frame
        p = tf.add_paragraph()
        p.text = "Results"
        p.font.size = Pt(28)
        left0 = 0.2
        top0 = 0.7
        height = Inches(2.5)
        height_w = Inches(1)
        #
        for i in range(0, 3):
            for k in range(0, 4):
                if imn >= num_ims:
                    break
                else:
                    left = Inches(left0 + k * 3.2)
                    top = Inches(top0 + i * 2.1)
                    slide.shapes.add_picture(img_names[imn], left, top, height=height)
                    txBox = slide.shapes.add_textbox(
                        left, top - Inches(0.2), width, height_w
                    )
                    tf = txBox.text_frame
                    p = tf.add_paragraph()

                    p.text = (
                        str(int(time[imn]))
                        + " s, "
                        + str("{:.2f}".format(coverage[imn]))
                        + " ML"
                    )
                    p.font.size = Pt(20)

                    imn = imn + 1

    # Third slide
    num_ims = len(hist_names)
    imn = 0
    while imn < num_ims:
        slide = prs.slides.add_slide(blank_slide_layout)
        width = height = Inches(1)
        top = Inches(-0.1)
        left = Inches(0.5)
        txBox = slide.shapes.add_textbox(left, top, width, height)
        tf = txBox.text_frame
        p = tf.add_paragraph()
        p.text = "Results: layer analysis"
        p.font.size = Pt(28)
        left0 = 0.2
        top0 = 0.7
        height = Inches(2)
        height_w = Inches(1)
        for i in range(0, 3):
            for k in range(0, 4):
                if imn >= num_ims:
                    break
                else:
                    left = Inches(left0 + k * 3.2)
                    top = Inches(top0 + i * 2.1)
                    slide.shapes.add_picture(
                        hist_names[imn], left, top + Inches(0.2), height=height
                    )

                    txBox = slide.shapes.add_textbox(
                        left, top - Inches(0.2), width, height_w
                    )
                    tf = txBox.text_frame

                    p = tf.add_paragraph()
                    p.text = (
                        str(int(time[imn]))
                        + " s, "
                        + str("{:.2f}".format(coverage[imn]))
                        + " ML"
                    )
                    p.font.size = Pt(20)

                    imn = imn + 1
    """
    slide = prs.slides.add_slide(blank_slide_layout)

    width = height = Inches(1)
    top = Inches(-0.1)
    left = Inches(0.5)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame

    p = tf.add_paragraph()
    p.text = "Results: Coverage dependence"
    p.font.size = Pt(28)

    height = Inches(5.5)
    top = Inches(1)
    left = Inches(0.7)

    file_name = entry_rec.get() + "_coverage" + ".png"

    slide.shapes.add_picture(file_name, left, top, height=height)
    """
    #
    slide = prs.slides.add_slide(blank_slide_layout)
    width = height = Inches(1)
    top = Inches(-0.1)
    left = Inches(0.5)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    p = tf.add_paragraph()
    p.text = "Number of events per deposition"
    p.font.size = Pt(28)
    height = Inches(5.5)
    top = Inches(1)
    left = Inches(0.7)
    slide.shapes.add_picture(
        dir_name + "Num_events_per_dep.png", left, top, height=height
    )
    #

    rec_num = 0
    while rec_num in (0, 1):
        try:
            prs.save(ppt_name)
            rec_num = 2
        except:
            if rec_num == 0:
                print("Close the ppt")
                rec_num = 1
            else:
                pass
