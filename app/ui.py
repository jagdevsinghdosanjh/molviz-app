import streamlit as st

def render_ui(geometry_data):
    st.subheader("📐 Molecular Geometry Details")
    st.write(f"**Geometry Type**: {geometry_data['geometry']}")
    st.write(f"**Bond Angles**: {geometry_data['angles']}°")
    st.write(f"**Bond Lengths**: {geometry_data['lengths']} Å")

    with st.expander("📊 Raw Data"):
        st.json(geometry_data)
